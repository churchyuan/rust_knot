use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::{KnotError, Result};
use crate::kmt::kmt_ring;
use crate::point::Point3;

/// Determine whether a 3_1 knot is a "two-arm" or "three-arm" state.
/// The logic mirrors `distinguish_2arm_and_3arm_states.m`:
/// 1. find exactly 3 projected crossings on the simplified open knot core,
/// 2. split the ring into 6 fragments using the 6 fractional bead indices,
/// 3. build polygons from the two fragments connecting the other two crossings,
/// 4. use MATLAB's ray-casting rule to classify the state.
pub fn get_31knot_arm_type(
    points: &[Point3],
    _core_points: &[Point3],
    _core_indices: &[usize],
    _table: &AlexanderTable,
    _config: &KnotConfig,
) -> Result<String> {
    if points.len() < 3 {
        return Err(KnotError::DataParse("Chain is too short".into()));
    }

    let mut last_error = None;
    for candidate in candidate_chains(points) {
        for projected in projected_variants(&candidate) {
            match classify_projected_points(&projected) {
                Ok(arm_type) => return Ok(arm_type),
                Err(err) => last_error = Some(err),
            }
        }
    }

    Err(last_error.unwrap_or_else(|| KnotError::DataParse("arm type classification failed".into())))
}

fn classify_projected_points(points: &[Point3]) -> Result<String> {
    let lc = points.len();
    if lc < 3 {
        return Err(KnotError::DataParse(
            "Chain is too short after simplification".into(),
        ));
    }
    let mut cross_xy: Vec<[f64; 2]> = Vec::new();
    let mut bead_index_of_eachcross: Vec<[f64; 2]> = Vec::new();

    for i in 0..lc {
        let a = points[i];
        let b = points[(i + 1) % lc];

        for j in (i + 2)..lc {
            if i == 0 && j + 1 == lc {
                continue;
            }

            let c = points[j];
            let d = points[(j + 1) % lc];

            let ax = b[0] - a[0];
            let bx = c[0] - d[0];
            let cx = c[0] - a[0];
            let dx = b[1] - a[1];
            let ex = c[1] - d[1];
            let fx = c[1] - a[1];
            let det = ax * ex - bx * dx;

            if det.abs() <= 1e-9 {
                continue;
            }

            let kab = (cx * ex - bx * fx) / det;
            let kcd = (ax * fx - cx * dx) / det;
            if !(kab > 0.0 && kab < 1.0 && kcd > 0.0 && kcd < 1.0) {
                continue;
            }

            let zab = a[2] * (1.0 - kab) + b[2] * kab;
            let zcd = c[2] * (1.0 - kcd) + d[2] * kcd;
            let ix = a[0] * (1.0 - kab) + b[0] * kab;
            let iy = a[1] * (1.0 - kab) + b[1] * kab;

            cross_xy.push([ix, iy]);

            if zab > zcd {
                bead_index_of_eachcross.push([i as f64 + kab, j as f64 + kcd]);
            } else {
                bead_index_of_eachcross.push([j as f64 + kcd, i as f64 + kab]);
            }
        }
    }

    if bead_index_of_eachcross.len() != 3 {
        return Err(KnotError::DataParse(format!(
            "Number of crossings is {}; not 3",
            bead_index_of_eachcross.len()
        )));
    }

    // PART 2: Sort three crossings and six fragment
    let mut bead_index = Vec::new();
    for row in &bead_index_of_eachcross {
        bead_index.push(row[0]);
        bead_index.push(row[1]);
    }
    bead_index.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut bead_index_of_eachfrag = vec![[0.0, 0.0]; 6];
    for i in 0..5 {
        bead_index_of_eachfrag[i] = [bead_index[i], bead_index[i + 1]];
    }
    bead_index_of_eachfrag[5] = [bead_index[5], bead_index[0] + lc as f64];

    let mut cross_index_of_eachcross = vec![0; 6];
    for i in 0..6 {
        let mut min_diff = f64::MAX;
        let mut best_cross = 0;
        for (c_idx, row) in bead_index_of_eachcross.iter().enumerate() {
            let d1 = (row[0] - bead_index[i]).abs();
            let d2 = (row[1] - bead_index[i]).abs();
            let d = d1.min(d2);
            if d < min_diff {
                min_diff = d;
                best_cross = c_idx + 1; // 1-based
            }
        }
        cross_index_of_eachcross[i] = best_cross;
    }

    let ci = &cross_index_of_eachcross;

    let mut cross_to_frag = [[[0; 2]; 3]; 3];
    let mut ntmp = [[0; 3]; 3];

    for i in 0..6 {
        let j = if i + 1 >= 6 { 0 } else { i + 1 };
        let c1 = ci[i] - 1;
        let c2 = ci[j] - 1;

        let m = ntmp[c1][c2];
        if m < 2 {
            cross_to_frag[c1][c2][m] = i + 1;
            ntmp[c1][c2] += 1;
        }
    }

    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                continue;
            }
            if cross_to_frag[i][j][0] == 0 {
                cross_to_frag[i][j][0] = cross_to_frag[j][i][0];
                cross_to_frag[i][j][1] = cross_to_frag[j][i][1];
            }
        }
    }

    let mut flag_twoarm = false;

    for icross in 1..=3 {
        let (ca, cb) = match icross {
            1 => (2, 3),
            2 => (3, 1),
            3 => (1, 2),
            _ => unreachable!(),
        };

        let frag_a = cross_to_frag[ca - 1][cb - 1][0];
        let frag_b = cross_to_frag[ca - 1][cb - 1][1];
        if frag_a == 0 || frag_b == 0 {
            return Err(KnotError::DataParse(format!(
                "Missing fragment mapping between crossings {} and {}",
                ca, cb
            )));
        }

        let frag_a_idx = frag_a - 1;
        let frag_b_idx = frag_b - 1;

        let ia = bead_index_of_eachfrag[frag_a_idx][0].ceil() as usize;
        let ib = bead_index_of_eachfrag[frag_a_idx][1].floor() as usize;
        let ic = bead_index_of_eachfrag[frag_b_idx][0].ceil() as usize;
        let id = bead_index_of_eachfrag[frag_b_idx][1].floor() as usize;

        let mut polygon = Vec::new();

        let mut k = ia;
        while k <= ib {
            let idx = k % lc;
            polygon.push([points[idx][0], points[idx][1]]);
            k += 1;
        }

        let mut k = id;
        while k >= ic {
            let idx = k % lc;
            polygon.push([points[idx][0], points[idx][1]]);
            if k == 0 {
                break;
            }
            k -= 1;
        }

        let lloop = polygon.len();
        if lloop < 3 {
            continue;
        }

        let mut nintersection = 0;
        let x0 = cross_xy[icross - 1][0];
        let y0 = cross_xy[icross - 1][1];

        for i in 0..lloop {
            let j = if i + 1 == lloop { 0 } else { i + 1 };
            let p1 = polygon[i];
            let p2 = polygon[j];

            if ((p1[1] - y0) * (p2[1] - y0)) < 0.0 && p2[0] > x0 {
                nintersection += 1;
            }
        }

        if nintersection % 2 == 1 {
            flag_twoarm = true;
            break;
        }
    }

    if flag_twoarm {
        Ok("two-arm".to_string())
    } else {
        Ok("three-arm".to_string())
    }
}

fn candidate_chains(points: &[Point3]) -> [Vec<Point3>; 2] {
    let mut simplified = points.to_vec();
    kmt_ring(&mut simplified);
    [points.to_vec(), simplified]
}

fn projected_variants(points: &[Point3]) -> [Vec<Point3>; 3] {
    [
        points.to_vec(),
        points.iter().copied().map(|p| [p[0], p[2], p[1]]).collect(),
        points.iter().copied().map(|p| [p[1], p[2], p[0]]).collect(),
    ]
}
