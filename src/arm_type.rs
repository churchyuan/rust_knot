use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::{KnotError, Result};
use crate::geometry::{cal_intersection, find_max_span, Intersection};
use crate::hull::hull_ends;
use crate::kmt::kmt_open_chain_with_indices;
use crate::knottype::get_knottype;
use crate::point::Point3;

/// Determine whether a 3_1 knot is a "two-arm" or "three-arm" state.
/// This algorithm should be run on the knot core (an open chain) of a 3_1 knot.
pub fn get_31knot_arm_type(
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String> {
    // 1. Verify that the input chain is indeed a 3_1 knot
    let knot_type = get_knottype(points, table, config)?;
    if knot_type != "3_1" {
        return Err(KnotError::DataParse(format!(
            "Expected a 3_1 knot, but found {}",
            knot_type
        )));
    }

    // 2. Initialize points with their original indices
    let mut working: Vec<(Point3, usize)> = points.iter().copied().enumerate().map(|(i, p)| (p, i)).collect();

    // 3. Apply KMT simplification while tracking original indices
    kmt_open_chain_with_indices(&mut working);

    // 4. Extend endpoints through convex hull to avoid closing segment interference
    if working.len() > 4 {
        let pts_only: Vec<Point3> = working.iter().map(|p| p.0).collect();
        if let Some((start_ext, end_ext)) =
            hull_ends(&pts_only, config.hull_plane_epsilon, config.extend_factor)
        {
            let first_idx = working.first().unwrap().1;
            let last_idx = working.last().unwrap().1;
            working.insert(0, (start_ext, first_idx));
            working.push((end_ext, last_idx));
        } else if config.debug {
            eprintln!("warning: hull_ends failed during arm type calculation, using infinity-closure fallback");
        }
    }

    let n = working.len();
    if n < 3 {
        return Err(KnotError::DataParse(
            "Simplified chain is too short to form crossings".into(),
        ));
    }

    // 5. Reorient to maximize 2D projection span
    let mut pts: Vec<Point3> = working.iter().map(|p| p.0).collect();
    find_max_span(&mut pts);

    // 6. Calculate all valid intersections
    let num_segs = n - 1;
    let mut intersections = vec![vec![Intersection::invalid(); num_segs]; num_segs];
    for i in 0..num_segs {
        for j in 0..num_segs {
            if i == j {
                continue;
            }
            intersections[i][j] = cal_intersection(&pts[i], &pts[i + 1], &pts[j], &pts[j + 1]);
        }
    }

    // Collect all unique crossings (segment pairs)
    let mut crossing_pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..num_segs {
        for j in (i + 1)..num_segs {
            if intersections[i][j].valid {
                crossing_pairs.push((i, j));
            }
        }
    }

    // We expect exactly 3 crossings for a minimal 3_1 knot projection after KMT
    if crossing_pairs.len() != 3 {
        return Err(KnotError::DataParse(format!(
            "Expected exactly 3 crossings for 3_1 knot, but found {}",
            crossing_pairs.len()
        )));
    }

    // 7. Map each crossing to its two original index locations on the chain
    // A segment `i` corresponds to the interval from `working[i].1` to `working[i+1].1`.
    // We approximate the crossing location by the starting index `working[i].1`.
    let mut crossings_locs = Vec::new();
    for (s1, s2) in crossing_pairs {
        let loc1 = working[s1].1;
        let loc2 = working[s2].1;
        // ensure loc1 < loc2 for consistency
        if loc1 < loc2 {
            crossings_locs.push([loc1, loc2]);
        } else {
            crossings_locs.push([loc2, loc1]);
        }
    }

    // 8. Order the 3 crossings by their first appearance on the chain
    // We collect all 6 locations and sort them to trace the chain traversal.
    let mut loc_to_crossing = Vec::new();
    for (c_idx, locs) in crossings_locs.iter().enumerate() {
        loc_to_crossing.push((locs[0], c_idx));
        loc_to_crossing.push((locs[1], c_idx));
    }
    loc_to_crossing.sort_by_key(|k| k.0);

    let mut ordered_crossings = Vec::new();
    for &(_, c_idx) in &loc_to_crossing {
        if !ordered_crossings.contains(&c_idx) {
            ordered_crossings.push(c_idx);
        }
    }

    if ordered_crossings.len() != 3 {
        return Err(KnotError::DataParse(
            "Failed to establish sequential order for 3 crossings".into(),
        ));
    }

    let c1_idx = ordered_crossings[0];
    let c2_idx = ordered_crossings[1];
    let c3_idx = ordered_crossings[2];

    let locs1 = crossings_locs[c1_idx];
    let locs2 = crossings_locs[c2_idx];
    let locs3 = crossings_locs[c3_idx];

    // 9. Evaluate the 4 paths between Crossing 1 and Crossing 3
    let paths = vec![
        (locs1[0], locs3[0]),
        (locs1[0], locs3[1]),
        (locs1[1], locs3[0]),
        (locs1[1], locs3[1]),
    ];

    // Calculate path lengths (physical distance based on original chain indices)
    let mut paths_with_dist: Vec<((usize, usize), usize)> = paths
        .into_iter()
        .map(|(a, b)| ((a, b), a.abs_diff(b)))
        .collect();

    // Sort by distance (shortest first)
    paths_with_dist.sort_by_key(|k| k.1);

    let shortest_1 = paths_with_dist[0].0;
    let shortest_2 = paths_with_dist[1].0;

    // Helper closure: checks if Crossing 2 passes through the given path
    // Meaning at least one of its two locations lies within the interval of the path.
    let passes_through = |path: (usize, usize), l2: [usize; 2]| -> bool {
        let min_p = path.0.min(path.1);
        let max_p = path.0.max(path.1);

        let l2a = l2[0];
        let l2b = l2[1];

        (l2a >= min_p && l2a <= max_p) || (l2b >= min_p && l2b <= max_p)
    };

    let pass1 = passes_through(shortest_1, locs2);
    let pass2 = passes_through(shortest_2, locs2);

    // 10. Decide arm type
    if pass1 && pass2 {
        Ok("two-arm".to_string())
    } else {
        Ok("three-arm".to_string())
    }
}
