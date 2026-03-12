use std::fs::File;
use std::io::BufReader;

use rust_knot::config::KnotConfig;
use rust_knot::hull::hull_ends;
use rust_knot::io::XyzFrameIter;
use rust_knot::kmt::kmt_open_chain;

fn run_file(path: &str, use_kmt: bool) {
    let file = File::open(path).expect("cannot open file");
    let reader = BufReader::new(file);
    let config = KnotConfig::default();

    let mut total = 0usize;
    let mut hull_ok = 0usize;
    let mut hull_fail = 0usize;
    let mut skipped_short = 0usize;

    for result in XyzFrameIter::new(reader) {
        let points = result.expect("parse error");
        total += 1;

        let mut working = points.clone();
        if use_kmt {
            kmt_open_chain(&mut working);
        }

        if working.len() <= 4 {
            skipped_short += 1;
            continue;
        }

        match hull_ends(&working, config.hull_plane_epsilon, config.extend_factor) {
            Some(_) => hull_ok += 1,
            None => {
                hull_fail += 1;
                if hull_fail <= 3 {
                    eprintln!(
                        "  frame {total}: hull_ends failed, chain len = {}",
                        working.len()
                    );
                }
            }
        }
    }

    let attempted = hull_ok + hull_fail;
    let pct = if attempted > 0 {
        hull_fail as f64 / attempted as f64 * 100.0
    } else {
        0.0
    };
    let kmt_label = if use_kmt { "WITH KMT" } else { "NO KMT" };
    println!(
        "{path:<45} [{kmt_label:<8}]  frames={total:<6} skip={skipped_short:<4} ok={hull_ok:<6} fail={hull_fail:<4} rate={pct:.2}%"
    );
}

fn main() {
    let files = [
        "test_manual/L300_knot3_1_open.xyz",
        "test_manual/L400_knot4_1_open.xyz",
    ];

    println!("=== Hull Fallback Statistics ===\n");
    for f in &files {
        run_file(f, true);
        run_file(f, false);
        println!();
    }
}
