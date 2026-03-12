use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::time::Instant;

use rust_knot::alexander_table::AlexanderTable;
use rust_knot::batch::process_frames_parallel;
use rust_knot::config::KnotConfig;
use rust_knot::io::read_data_xyz_frames;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <table_file> <xyz_file> [target_type] [--ring] [--fast] [--debug] [--output <path>]",
            args[0]
        );
        eprintln!("  table_file:    path to table_knot_Alexander_polynomial.txt");
        eprintln!("  xyz_file:      path to XYZ coordinate file (single or multi-frame)");
        eprintln!("  target_type:   (optional) knot type to search for core, e.g. '3_1'");
        eprintln!("  --ring:        treat chain as a closed ring");
        eprintln!("  --fast:        enable KMT simplification");
        eprintln!("  --debug:       enable debug output");
        eprintln!("  --output <path>: write knot_index log (default: knot_index.txt)");
        std::process::exit(1);
    }

    let table_path = &args[1];
    let xyz_path = &args[2];

    // Parse flags
    let mut config = KnotConfig {
        faster: true,
        ..KnotConfig::default()
    };
    let mut target_type: Option<String> = None;
    let mut output_path = String::from("knot_index.txt");
    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--ring" => config.is_ring = true,
            "--fast" => config.faster = true,
            "--debug" => config.debug = true,
            "--output" => {
                i += 1;
                if i < args.len() {
                    output_path = args[i].clone();
                }
            }
            s if !s.starts_with("--") => target_type = Some(s.to_string()),
            other => eprintln!("warning: unknown flag '{other}'"),
        }
        i += 1;
    }

    // Load Alexander polynomial table
    let t0 = Instant::now();
    let table = AlexanderTable::from_file(table_path).expect("failed to load Alexander table");
    eprintln!(
        "Loaded Alexander table ({} entries) in {:?}",
        table.len(),
        t0.elapsed()
    );

    // Read all frames
    let t0 = Instant::now();
    let file = File::open(xyz_path).expect("failed to open XYZ file");
    let mut reader = BufReader::new(file);
    let frames = read_data_xyz_frames(&mut reader).expect("failed to parse XYZ file");

    if frames.is_empty() {
        eprintln!("No frames found in {xyz_path}");
        std::process::exit(1);
    }

    eprintln!(
        "Read {} frame(s) ({} points/frame) in {:?}",
        frames.len(),
        frames[0].len(),
        t0.elapsed()
    );
    eprintln!(
        "Config: is_ring={}, faster={}, extend_factor={}, num_rotations={}",
        config.is_ring, config.faster, config.extend_factor, config.num_rotations
    );

    // Process all frames in parallel
    let t0 = Instant::now();
    let results =
        process_frames_parallel(&frames, &table, &config, target_type.as_deref());
    let compute_time = t0.elapsed();

    // Write knot_index.txt (sequential, in frame order)
    let out_file = File::create(&output_path).expect("failed to create output file");
    let mut writer = BufWriter::new(out_file);

    writeln!(writer, "# frame\tknottype\tknot_start\tknot_end\tknot_size")
        .expect("write header failed");

    for r in &results {
        if let Some(ref err) = r.error {
            eprintln!("frame {}: error: {err}", r.frame);
        }
        let ktype = if r.error.is_some() { "ERROR" } else { &r.knot_type };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            r.frame, ktype, r.knot_start, r.knot_end, r.knot_size
        )
        .expect("write failed");
    }
    writer.flush().expect("flush failed");

    // Summary to stderr
    let n_knotted = results.iter().filter(|r| r.knot_type != "1" && r.error.is_none()).count();
    let n_errors = results.iter().filter(|r| r.error.is_some()).count();
    eprintln!(
        "Done: {} frames, {} knotted, {} unknot, {} errors — {:?} compute, written to {}",
        results.len(),
        n_knotted,
        results.len() - n_knotted - n_errors,
        n_errors,
        compute_time,
        output_path
    );

    // Single frame: also print details to stdout for convenience
    if results.len() == 1 {
        let r = &results[0];
        if let Some(ref err) = r.error {
            println!("Error: {err}");
        } else {
            println!("Knot type: {}", r.knot_type);
            if r.knot_start >= 0 {
                println!(
                    "Knot core: [{}, {}], size = {}",
                    r.knot_start, r.knot_end, r.knot_size
                );
            }
        }
    }
}
