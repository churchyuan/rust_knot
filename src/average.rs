use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::{KnotError, Result};
use crate::io::write_data_xyz;
use crate::point::Point3;

#[derive(Clone, Debug)]
pub struct AlignResult {
    pub aligned: Vec<Point3>,
    pub rmsd: f64,
    pub shift: usize,
    pub reversed: bool,
}

fn mean_point(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let mut m = [0.0_f64; 3];
    for p in points {
        m[0] += p[0];
        m[1] += p[1];
        m[2] += p[2];
    }
    [m[0] / n, m[1] / n, m[2] / n]
}

fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: Point3, b: Point3) -> Point3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn mat_vec(m: [[f64; 3]; 3], v: Point3) -> Point3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn quat_normalize(mut q: [f64; 4]) -> [f64; 4] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n == 0.0 {
        q[0] = 1.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 0.0;
        return q;
    }
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}

fn quat_to_rot(q: [f64; 4]) -> [[f64; 3]; 3] {
    let w = q[0];
    let x = q[1];
    let y = q[2];
    let z = q[3];
    let ww = w * w;
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    [
        [ww + xx - yy - zz, 2.0 * (xy - wz), 2.0 * (xz + wy)],
        [2.0 * (xy + wz), ww - xx + yy - zz, 2.0 * (yz - wx)],
        [2.0 * (xz - wy), 2.0 * (yz + wx), ww - xx - yy + zz],
    ]
}

fn kabsch_align(mobile: &[Point3], target: &[Point3]) -> Result<(Vec<Point3>, f64)> {
    if mobile.len() != target.len() {
        return Err(KnotError::DataParse(format!(
            "point count mismatch: mobile {}, target {}",
            mobile.len(),
            target.len()
        )));
    }
    if mobile.is_empty() {
        return Err(KnotError::EmptyChain);
    }

    let m0 = mean_point(mobile);
    let t0 = mean_point(target);

    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut sxz = 0.0;
    let mut syx = 0.0;
    let mut syy = 0.0;
    let mut syz = 0.0;
    let mut szx = 0.0;
    let mut szy = 0.0;
    let mut szz = 0.0;

    for (pm, pt) in mobile.iter().zip(target.iter()) {
        let p = sub(*pm, m0);
        let q = sub(*pt, t0);
        sxx += p[0] * q[0];
        sxy += p[0] * q[1];
        sxz += p[0] * q[2];
        syx += p[1] * q[0];
        syy += p[1] * q[1];
        syz += p[1] * q[2];
        szx += p[2] * q[0];
        szy += p[2] * q[1];
        szz += p[2] * q[2];
    }

    let tr = sxx + syy + szz;

    let n = [
        [tr, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];

    let mut q = [1.0_f64, 0.0, 0.0, 0.0];
    for _ in 0..64 {
        let qn = [
            n[0][0] * q[0] + n[0][1] * q[1] + n[0][2] * q[2] + n[0][3] * q[3],
            n[1][0] * q[0] + n[1][1] * q[1] + n[1][2] * q[2] + n[1][3] * q[3],
            n[2][0] * q[0] + n[2][1] * q[1] + n[2][2] * q[2] + n[2][3] * q[3],
            n[3][0] * q[0] + n[3][1] * q[1] + n[3][2] * q[2] + n[3][3] * q[3],
        ];
        q = quat_normalize(qn);
    }

    let r = quat_to_rot(q);

    let mut aligned = Vec::with_capacity(mobile.len());
    let mut sum_sq = 0.0_f64;
    for (pm, pt) in mobile.iter().zip(target.iter()) {
        let p = sub(*pm, m0);
        let pr = mat_vec(r, p);
        let a = add(pr, t0);
        let d = sub(a, *pt);
        sum_sq += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        aligned.push(a);
    }
    let rmsd = (sum_sq / mobile.len() as f64).sqrt();
    Ok((aligned, rmsd))
}

fn shift_points(points: &[Point3], shift: usize) -> Vec<Point3> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(points[(i + shift) % n]);
    }
    out
}

fn best_cyclic_align(mobile: &[Point3], target: &[Point3], allow_reverse: bool) -> Result<AlignResult> {
    let n = mobile.len();
    if n != target.len() {
        return Err(KnotError::DataParse(format!(
            "point count mismatch: mobile {}, target {}",
            n,
            target.len()
        )));
    }
    if n == 0 {
        return Err(KnotError::EmptyChain);
    }

    let mut best: Option<AlignResult> = None;
    for k in 0..n {
        let shifted = shift_points(mobile, k);
        let (aligned, rmsd) = kabsch_align(&shifted, target)?;
        let cand = AlignResult {
            aligned,
            rmsd,
            shift: k,
            reversed: false,
        };
        if best.as_ref().map(|b| cand.rmsd < b.rmsd).unwrap_or(true) {
            best = Some(cand);
        }
    }

    if allow_reverse {
        let mut rev = mobile.to_vec();
        rev.reverse();
        for k in 0..n {
            let shifted = shift_points(&rev, k);
            let (aligned, rmsd) = kabsch_align(&shifted, target)?;
            let cand = AlignResult {
                aligned,
                rmsd,
                shift: k,
                reversed: true,
            };
            if best.as_ref().map(|b| cand.rmsd < b.rmsd).unwrap_or(true) {
                best = Some(cand);
            }
        }
    }

    best.ok_or_else(|| KnotError::DataParse("failed to align".into()))
}

fn align(mobile: &[Point3], target: &[Point3], cyclic: bool, allow_reverse: bool) -> Result<AlignResult> {
    if cyclic {
        return best_cyclic_align(mobile, target, allow_reverse);
    }
    let (aligned, rmsd) = kabsch_align(mobile, target)?;
    Ok(AlignResult {
        aligned,
        rmsd,
        shift: 0,
        reversed: false,
    })
}

pub fn compute_average(
    frames: &[Vec<Point3>],
    reference: &[Point3],
    cyclic: bool,
    allow_reverse: bool,
    max_iter: usize,
    tol: f64,
) -> Result<(Vec<Point3>, Vec<AlignResult>)> {
    if frames.is_empty() {
        return Err(KnotError::DataParse("no frames provided".into()));
    }
    let n = reference.len();
    if n == 0 {
        return Err(KnotError::EmptyChain);
    }
    for (i, f) in frames.iter().enumerate() {
        if f.len() != n {
            return Err(KnotError::DataParse(format!(
                "frame {} point count mismatch: expected {}, got {}",
                i,
                n,
                f.len()
            )));
        }
    }

    let mut avg = reference.to_vec();
    let max_iter = max_iter.max(1);
    let tol = tol.max(0.0);

    for _ in 0..max_iter {
        let mut aligned_all: Vec<Vec<Point3>> = Vec::with_capacity(frames.len());
        for fr in frames {
            let r = align(fr, &avg, cyclic, allow_reverse)?;
            aligned_all.push(r.aligned);
        }

        let mut new_avg = vec![[0.0_f64; 3]; n];
        let inv = 1.0 / frames.len() as f64;
        for a in &aligned_all {
            for (i, p) in a.iter().enumerate() {
                new_avg[i][0] += p[0] * inv;
                new_avg[i][1] += p[1] * inv;
                new_avg[i][2] += p[2] * inv;
            }
        }

        let mut sum_sq = 0.0_f64;
        for (a, b) in new_avg.iter().zip(avg.iter()) {
            let d = sub(*a, *b);
            sum_sq += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
        let delta = (sum_sq / n as f64).sqrt();
        avg = new_avg;
        if delta < tol {
            break;
        }
    }

    let mut final_results = Vec::with_capacity(frames.len());
    for fr in frames {
        final_results.push(align(fr, &avg, cyclic, allow_reverse)?);
    }
    Ok((avg, final_results))
}

pub fn write_average_outputs(
    outdir: &Path,
    avg: &[Point3],
    aligned_results: &[AlignResult],
    selected_frames: &[usize],
) -> Result<()> {
    std::fs::create_dir_all(outdir)?;

    let mut avg_file = BufWriter::new(File::create(outdir.join("average.xyz"))?);
    write_data_xyz(avg, &mut avg_file)?;
    avg_file.flush().map_err(KnotError::Io)?;

    let mut sum_rmsd2 = 0.0_f64;
    for r in aligned_results {
        sum_rmsd2 += r.rmsd * r.rmsd;
    }
    let tube_rmsd_mean = if aligned_results.is_empty() {
        f64::NAN
    } else {
        (sum_rmsd2 / aligned_results.len() as f64).sqrt()
    };

    let n = avg.len();
    let mut msd_atom = vec![0.0_f64; n];
    if !aligned_results.is_empty() {
        let inv = 1.0 / aligned_results.len() as f64;
        for r in aligned_results {
            for i in 0..n {
                let d = sub(r.aligned[i], avg[i]);
                msd_atom[i] += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) * inv;
            }
        }
    }
    let mut msd_sum = 0.0_f64;
    let mut msd_max = 0.0_f64;
    for v in &msd_atom {
        msd_sum += *v;
        if *v > msd_max {
            msd_max = *v;
        }
    }
    let tube_atom_rms = if n == 0 { f64::NAN } else { (msd_sum / n as f64).sqrt() };
    let tube_atom_max = msd_max.sqrt();

    let mut summary = BufWriter::new(File::create(outdir.join("summary.txt"))?);
    writeln!(summary, "frames {}", aligned_results.len()).map_err(KnotError::Io)?;
    writeln!(summary, "points {}", avg.len()).map_err(KnotError::Io)?;
    writeln!(summary, "tube_rmsd_mean {:.10e}", tube_rmsd_mean).map_err(KnotError::Io)?;
    writeln!(summary, "tube_atom_rms {:.10e}", tube_atom_rms).map_err(KnotError::Io)?;
    writeln!(summary, "tube_atom_max {:.10e}", tube_atom_max).map_err(KnotError::Io)?;
    summary.flush().map_err(KnotError::Io)?;

    let mut sel = BufWriter::new(File::create(outdir.join("selected_frames.csv"))?);
    writeln!(sel, "index,frame").map_err(KnotError::Io)?;
    for (i, fr) in selected_frames.iter().enumerate() {
        writeln!(sel, "{},{}", i, fr).map_err(KnotError::Io)?;
    }
    sel.flush().map_err(KnotError::Io)?;

    let mut rmsd = BufWriter::new(File::create(outdir.join("rmsd.csv"))?);
    writeln!(rmsd, "index,frame,rmsd,shift,reversed").map_err(KnotError::Io)?;
    for (i, (fr, r)) in selected_frames.iter().zip(aligned_results.iter()).enumerate() {
        writeln!(
            rmsd,
            "{},{},{:.10e},{},{}",
            i,
            fr,
            r.rmsd,
            r.shift,
            if r.reversed { 1 } else { 0 }
        )
        .map_err(KnotError::Io)?;
    }
    rmsd.flush().map_err(KnotError::Io)?;

    let mut msd = BufWriter::new(File::create(outdir.join("msd_atom.csv"))?);
    writeln!(msd, "atom_index,msd,rms").map_err(KnotError::Io)?;
    for (i, v) in msd_atom.iter().enumerate() {
        writeln!(msd, "{},{:.10e},{:.10e}", i, v, v.sqrt()).map_err(KnotError::Io)?;
    }
    msd.flush().map_err(KnotError::Io)?;

    Ok(())
}
