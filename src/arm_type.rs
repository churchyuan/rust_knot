use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::{KnotError, Result};
use crate::geometry::find_max_span;
use crate::point::Point3;

use crate::kmt::kmt_open_chain_with_indices;

/// Determine whether a 3_1 knot is a "two-arm" or "three-arm" state.
/// This algorithm should be run on the FULL chain of a 3_1 knot,
/// to correctly handle ring wrap-around and precise geometry intersections.
pub fn get_31knot_arm_type(
    points: &[Point3],
    core_points: &[Point3],
    core_indices: &[usize],
    _table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String> {
    let lc = points.len();
    if lc < 3 {
        return Err(KnotError::DataParse("Chain is too short".into()));
    }

    let mut working_kmt: Vec<(Point3, usize)> = core_points.iter().copied().zip(core_indices.iter().copied()).collect();
    kmt_open_chain_with_indices(&mut working_kmt);

    use crate::hull::hull_ends;
    if working_kmt.len() > 4 {
        let pts_only: Vec<Point3> = working_kmt.iter().map(|p| p.0).collect();
        if let Some((start_ext, end_ext)) =
            hull_ends(&pts_only, config.hull_plane_epsilon, config.extend_factor)
        {
            let first_idx = working_kmt.first().unwrap().1;
            let last_idx = working_kmt.last().unwrap().1;
            working_kmt.insert(0, (start_ext, first_idx));
            working_kmt.push((end_ext, last_idx));
        }
    }
    let mut working_pts: Vec<Point3> = working_kmt.iter().map(|p| p.0).collect();
    find_max_span(&mut working_pts);

    // Apply the exact same rotation to full points
    let mut working_pts_full = points.to_vec();
    find_max_span(&mut working_pts_full);

    // If the input chain has more than 3 crossings, the MATLAB logic won't work.
    // The previous simple logic worked beautifully for open chains even when they had >3 crossings?
    // Wait, in previous commit, it found EXACTLY 3 crossings!
    // Why did it find 3 crossings before? Let's trace it carefully.
    // 1. `kmt_open_chain_with_indices` on `working: Vec<(Point3, usize)> = points...`
    // 2. `hull_ends` and add `start_ext`, `end_ext` to `working`.
    // 3. `let mut working_pts: Vec<Point3> = working.iter().map(|p| p.0).collect();`
    // 4. `find_max_span(&mut working_pts);`
    // 5. Loop `num_segs = n - 1;` where `n = working_pts.len()`.
    // 6. Double loop `i = 0..num_segs` and `j = 0..num_segs`.
    // 7. `cal_intersection` and filter `valid`.
    // 8. Output was 3 crossings exactly!
    //
    // NOW what are we doing differently?
    // 1. `let mut working_kmt: Vec<(Point3, usize)> = points...`
    // 2. `kmt_open_chain_with_indices(&mut working_kmt);`
    // 3. `hull_ends`
    // 4. `let mut working_pts: Vec<Point3> = working_kmt.iter().map(|p| p.0).collect();`
    // 5. `find_max_span(&mut working_pts);`
    // 6. Loop `num_segs = n - 1;` where `n = working_pts.len()`.
    // 7. `for i in 0..num_segs`, `for j in (i + 2)..num_segs`.
    // 8. Output is 8 crossings!
    // Wait, `for j in 0..num_segs` with `i < j` is the same as `for j in (i + 1)..num_segs`.
    // But we are using `cal_intersection` on `&a, &b, &c, &d` where a=i, b=i+1, c=j, d=j+1.
    // Wait, if `for j in 0..num_segs` with `i < j` gave 3 crossings before...
    // Let me check my previous commit AGAIN.
    // Did it find 3 crossings? YES, because `test_31knot_arm_type` passed!
    // Wait, maybe it panicked but the test was marked `#[ignore]`? No, it passed!
    // Let's use `crate::geometry::cal_intersection` exactly like before.

    let n = working_pts.len();
    let num_segs = n - 1;
    let mut intersections = vec![vec![crate::geometry::Intersection::invalid(); num_segs]; num_segs];
    for i in 0..num_segs {
        for j in 0..num_segs {
            if i == j {
                continue;
            }
            intersections[i][j] = crate::geometry::cal_intersection(
                &working_pts[i],
                &working_pts[i + 1],
                &working_pts[j],
                &working_pts[j + 1],
            );
        }
    }

    let mut crossing_pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..num_segs {
        for j in (i + 1)..num_segs {
            if intersections[i][j].valid {
                crossing_pairs.push((i, j));
            }
        }
    }

    let ncross = crossing_pairs.len();
    if ncross != 3 {
        return Err(KnotError::DataParse(format!(
            "Number of crossings is {}; not 3",
            ncross
        )));
    }

    let mut cross_xy: Vec<[f64; 2]> = Vec::new();
    let mut bead_index_of_eachcross: Vec<[f64; 2]> = Vec::new();

    for (i, j) in crossing_pairs {
        let a = working_pts[i];
        let b = working_pts[i + 1];
        let c = working_pts[j];
        let d = working_pts[j + 1];

        let ax = b[0] - a[0];
        let bx = c[0] - d[0];
        let cx = c[0] - a[0];
        let dx = b[1] - a[1];
        let ex = c[1] - d[1];
        let fx = c[1] - a[1];
        let det = ax * ex - bx * dx;
        
        let mut kab = 0.5;
        let mut kcd = 0.5;
        if det.abs() > 1e-9 {
            kab = (cx * ex - bx * fx) / det;
            kcd = (ax * fx - cx * dx) / det;
        }

        let zab = a[2] * (1.0 - kab) + b[2] * kab;
        let zcd = c[2] * (1.0 - kcd) + d[2] * kcd;
        
        let ix = a[0] * (1.0 - kab) + b[0] * kab;
        let iy = a[1] * (1.0 - kab) + b[1] * kab;
        cross_xy.push([ix, iy]);

        let orig_i = working_kmt[i].1 as f64;
        let orig_j = working_kmt[j].1 as f64;

        if zab > zcd {
            bead_index_of_eachcross.push([orig_i + kab, orig_j + kcd]);
        } else {
            bead_index_of_eachcross.push([orig_j + kcd, orig_i + kab]);
        }
    }

    // PART 2: Sort three crossings and six fragment
    let mut bead_index = Vec::new();
    for row in &bead_index_of_eachcross {
        bead_index.push(row[0]);
        bead_index.push(row[1]);
    }
    // Sort all 6 fractional indices
    bead_index.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut bead_index_of_eachfrag = vec![[0.0, 0.0]; 6];
    for i in 0..5 {
        bead_index_of_eachfrag[i] = [bead_index[i], bead_index[i + 1]];
    }
    bead_index_of_eachfrag[5] = [bead_index[5], bead_index[0] + lc as f64];

    // Map each sorted location back to its crossing (1, 2, or 3) -> 1-based index to match MATLAB
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

    // CrossIndex_to_FragIndex: Maps pair of crossing indices (e.g. 1 and 3) to the two Fragments that connect them.
    // MATLAB array is 3x3x2. Rust: [[[0; 2]; 3]; 3]
    let mut cross_to_frag = [[[0; 2]; 3]; 3];
    let mut ntmp = [[0; 3]; 3];

    for i in 0..6 {
        let j = if i + 1 >= 6 { 0 } else { i + 1 };
        let c1 = ci[i] - 1; // 0-based
        let c2 = ci[j] - 1; // 0-based
        
        let m = ntmp[c1][c2];
        if m < 2 {
            cross_to_frag[c1][c2][m] = i + 1; // 1-based fragment index
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

    // PART 3: Ray casting algorithm for PIP analysis
    let mut flag_twoarm = false;

    // Now we must ensure that `points` is projected to the exact SAME XY plane
    // so that when we extract `polygon` from `points`, the coordinates match `cross_xy`.
    // Wait, `cross_xy` was calculated from `working_pts` (which had `find_max_span` applied).
    // So we must apply `find_max_span` to `points` BEFORE we extract the polygon!
    let mut working_pts_full = points.to_vec();
    find_max_span(&mut working_pts_full);

    for icross in 1..=3 {
        let (ca, cb) = match icross {
            1 => (2, 3),
            2 => (3, 1),
            3 => (1, 2),
            _ => unreachable!(),
        };

        // 1-based indices from MATLAB logic -> 0-based for array access
        let frag_a_idx = cross_to_frag[ca - 1][cb - 1][0] - 1;
        let frag_b_idx = cross_to_frag[ca - 1][cb - 1][1] - 1;

        let ia = bead_index_of_eachfrag[frag_a_idx][0].ceil() as usize;
        let ib = bead_index_of_eachfrag[frag_a_idx][1].floor() as usize;
        let ic = bead_index_of_eachfrag[frag_b_idx][0].ceil() as usize;
        let id = bead_index_of_eachfrag[frag_b_idx][1].floor() as usize;

        let mut polygon = Vec::new();

        // Forward loop for first fragment
        let mut k = ia;
        while k <= ib {
            let idx = if k >= lc { k - lc } else { k };
            polygon.push([working_pts_full[idx][0], working_pts_full[idx][1]]);
            k += 1;
        }

        // Backward loop for second fragment
        let mut k = id;
        while k >= ic {
            let idx = if k >= lc { k - lc } else { k };
            polygon.push([working_pts_full[idx][0], working_pts_full[idx][1]]);
            if k == 0 {
                break;
            }
            k -= 1;
        }

        let lloop = polygon.len();
        if lloop < 3 {
            continue; // Not enough points for a polygon
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
