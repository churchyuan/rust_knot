use std::fs::File;
use std::io::BufReader;

use rust_knot::batch::process_frame;
use rust_knot::config::KnotConfig;
use rust_knot::io::LammpsFrameIter;
use rust_knot::AlexanderTable;

fn evaluate_dataset(path: &str, expected_arm_type: &str) -> (usize, usize, Vec<String>) {
    let file = File::open(path).unwrap_or_else(|e| panic!("failed to open {path}: {e}"));
    let reader = BufReader::new(file);
    let table = AlexanderTable::builtin();
    let config = KnotConfig {
        is_ring: true,
        faster: true,
        ..KnotConfig::default()
    };

    let mut total = 0usize;
    let mut matched = 0usize;
    let mut failures = Vec::new();

    for (frame_idx, frame_result) in LammpsFrameIter::new(reader).enumerate() {
        let points =
            frame_result.unwrap_or_else(|e| panic!("failed to read frame {frame_idx}: {e}"));
        let result = process_frame(frame_idx, &points, &table, &config, Some("3_1"), true);
        total += 1;

        if result.error.is_none()
            && result.knot_type == "3_1"
            && result.arm_type.as_deref() == Some(expected_arm_type)
        {
            matched += 1;
        } else {
            failures.push(format!(
                "frame {} => knot_type={}, arm_type={:?}, error={:?}, warnings={:?}",
                frame_idx, result.knot_type, result.arm_type, result.error, result.warnings
            ));
        }
    }

    (matched, total, failures)
}

#[test]
fn test_two_arm_lammps_accuracy() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/twothree_arm_test/ALL_two_arm.lammpstrj"
    );
    let (matched, total, failures) = evaluate_dataset(path, "two-arm");
    eprintln!("two-arm accuracy: {matched}/{total}");
    for failure in failures.iter().take(10) {
        eprintln!("{failure}");
    }
    assert!(
        matched + 1 >= total,
        "two-arm dataset accuracy regressed: {matched}/{total}"
    );
}

#[test]
fn test_three_arm_lammps_accuracy() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/twothree_arm_test/ALL_three_arm.lammpstrj"
    );
    let (matched, total, failures) = evaluate_dataset(path, "three-arm");
    eprintln!("three-arm accuracy: {matched}/{total}");
    for failure in failures.iter().take(10) {
        eprintln!("{failure}");
    }
    assert!(
        matched == total,
        "three-arm dataset should classify every frame: {matched}/{total}"
    );
}
