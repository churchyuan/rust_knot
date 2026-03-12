use rayon::prelude::*;

use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::knotsize::find_knot_core;
use crate::knottype::get_knottype;
use crate::point::Point3;

/// Result of processing a single frame.
#[derive(Clone, Debug)]
pub struct FrameResult {
    pub frame: usize,
    pub knot_type: String,
    pub knot_start: i32,
    pub knot_end: i32,
    pub knot_size: i32,
    pub error: Option<String>,
}

/// Process a single frame: identify knot type, then find knot core.
///
/// If `target_type` is provided, the core search uses it instead of the
/// detected knot type.
pub fn process_frame(
    frame: usize,
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
) -> FrameResult {
    let knot_type = match get_knottype(points, table, config) {
        Ok(t) => t,
        Err(e) => {
            return FrameResult {
                frame,
                knot_type: String::new(),
                knot_start: -1,
                knot_end: -1,
                knot_size: 0,
                error: Some(e.to_string()),
            };
        }
    };

    let (start, end, size) = if knot_type == "1" {
        (-1, -1, 0)
    } else {
        let search = target_type.unwrap_or(&knot_type);
        match find_knot_core(points, search, table, config) {
            Ok(core) if core.matched => (core.left, core.right, core.size),
            _ => (-1, -1, 0),
        }
    };

    FrameResult {
        frame,
        knot_type,
        knot_start: start,
        knot_end: end,
        knot_size: size,
        error: None,
    }
}

/// Process all frames in parallel using rayon, returning results in frame order.
pub fn process_frames_parallel(
    frames: &[Vec<Point3>],
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
) -> Vec<FrameResult> {
    frames
        .par_iter()
        .enumerate()
        .map(|(i, points)| process_frame(i, points, table, config, target_type))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point3;

    fn straight_line(n: usize) -> Vec<Point3> {
        (0..n).map(|i| [i as f64, 0.0, 0.0]).collect()
    }

    #[test]
    fn test_process_frame_unknot() {
        let table_data = "0_1\t1\n3_1\t-1+1*t\n";
        let table =
            AlexanderTable::from_reader(std::io::Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let pts = straight_line(50);
        let result = process_frame(0, &pts, &table, &config, None);
        assert_eq!(result.knot_type, "1");
        assert_eq!(result.knot_start, -1);
        assert!(result.error.is_none());
    }

    #[test]
    fn test_process_frames_parallel_ordering() {
        let table_data = "0_1\t1\n";
        let table =
            AlexanderTable::from_reader(std::io::Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let frames: Vec<Vec<Point3>> = (0..5).map(|_| straight_line(50)).collect();
        let results = process_frames_parallel(&frames, &table, &config, None);
        assert_eq!(results.len(), 5);
        for (i, r) in results.iter().enumerate() {
            assert_eq!(r.frame, i);
        }
    }
}
