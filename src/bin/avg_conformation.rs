use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use rust_knot::average::{compute_average, write_average_outputs};
use rust_knot::error::{KnotError, Result};
use rust_knot::io::{read_data_xyz, LammpsFrameIter, XyzFrameIter};
use rust_knot::point::Point3;
use rust_knot::selection::{index_by_frame, read_knot_index, SelectionRow};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputFormat {
    Xyz,
    Lammps,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Mode {
    FullChain,
    CoreByLength,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ArmState {
    TwoArm,
    ThreeArm,
    Both,
}

fn print_usage(prog: &str) {
    eprintln!(
        "Usage: {prog} <input_file> --selection <knot_index.txt> [--format <xyz|lammps>] [--mode <full-chain|core-by-length>] [--target-type <knot>] [--arm-state <two-arm|three-arm|both>] [--reference <path>] [--outdir <path>] [--no-cyclic] [--no-reverse] [--max-iter <n>] [--tol <x>]"
    );
}

fn require_arg(args: &[String], i: usize, flag: &str) -> String {
    if i >= args.len() {
        eprintln!("error: {flag} requires a value");
        std::process::exit(1);
    }
    args[i].clone()
}

fn parse_input_format(value: &str) -> InputFormat {
    match value.to_ascii_lowercase().as_str() {
        "xyz" => InputFormat::Xyz,
        "lammps" | "lmp" | "lammpstrj" | "dump" => InputFormat::Lammps,
        _ => {
            eprintln!("error: unsupported --format '{value}', expected 'xyz' or 'lammps'");
            std::process::exit(1);
        }
    }
}

fn detect_input_format(path: &str) -> InputFormat {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".lammpstrj") || lower.ends_with(".lmp") || lower.ends_with(".dump") {
        InputFormat::Lammps
    } else {
        InputFormat::Xyz
    }
}

fn parse_mode(value: &str) -> Mode {
    match value.to_ascii_lowercase().as_str() {
        "full-chain" | "full" => Mode::FullChain,
        "core-by-length" | "core" => Mode::CoreByLength,
        _ => {
            eprintln!(
                "error: unsupported --mode '{value}', expected 'full-chain' or 'core-by-length'"
            );
            std::process::exit(1);
        }
    }
}

fn parse_arm_state(value: &str) -> ArmState {
    match value.to_ascii_lowercase().as_str() {
        "two-arm" | "two" => ArmState::TwoArm,
        "three-arm" | "three" => ArmState::ThreeArm,
        "both" => ArmState::Both,
        _ => {
            eprintln!("error: unsupported --arm-state '{value}', expected 'two-arm', 'three-arm', or 'both'");
            std::process::exit(1);
        }
    }
}

fn load_reference(path: &Path) -> Result<Vec<Point3>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    read_data_xyz(&mut reader)
}

fn extract_core(points: &[Point3], row: &SelectionRow) -> Result<Vec<Point3>> {
    if row.knot_start < 0 || row.knot_end < 0 {
        return Err(KnotError::DataParse(format!(
            "frame {} has no knot core indices",
            row.frame
        )));
    }
    let n = points.len() as i32;
    if n <= 0 {
        return Err(KnotError::EmptyChain);
    }
    let mut curr = row.knot_start.rem_euclid(n);
    let right = row.knot_end.rem_euclid(n);
    let mut out = Vec::new();
    loop {
        out.push(points[curr as usize]);
        if curr == right {
            break;
        }
        curr = (curr + 1) % n;
    }
    Ok(out)
}

fn matches_arm(row: &SelectionRow, arm_state: ArmState) -> Result<bool> {
    match arm_state {
        ArmState::Both => Ok(true),
        ArmState::TwoArm => Ok(row
            .arm_type
            .as_deref()
            .ok_or_else(|| KnotError::DataParse("selection file missing arm_type column".into()))?
            == "two-arm"),
        ArmState::ThreeArm => Ok(row
            .arm_type
            .as_deref()
            .ok_or_else(|| KnotError::DataParse("selection file missing arm_type column".into()))?
            == "three-arm"),
    }
}

fn run_cli(args: &[String], prog: &str) -> Result<()> {
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage(prog);
        std::process::exit(if args.is_empty() { 1 } else { 0 });
    }

    let input_path = &args[0];
    let mut input_format: Option<InputFormat> = None;
    let mut selection_path: Option<PathBuf> = None;
    let mut mode = Mode::FullChain;
    let mut target_type = "3_1".to_string();
    let mut arm_state = ArmState::Both;
    let mut reference_path: Option<PathBuf> = None;
    let mut outdir = PathBuf::from("avg_out");
    let mut no_cyclic = false;
    let mut no_reverse = false;
    let mut max_iter = 10usize;
    let mut tol = 1e-8_f64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--format" => {
                i += 1;
                let v = require_arg(args, i, "--format");
                input_format = Some(parse_input_format(&v));
            }
            "--selection" => {
                i += 1;
                selection_path = Some(PathBuf::from(require_arg(args, i, "--selection")));
            }
            "--mode" => {
                i += 1;
                let v = require_arg(args, i, "--mode");
                mode = parse_mode(&v);
            }
            "--target-type" => {
                i += 1;
                target_type = require_arg(args, i, "--target-type");
            }
            "--arm-state" => {
                i += 1;
                let v = require_arg(args, i, "--arm-state");
                arm_state = parse_arm_state(&v);
            }
            "--reference" => {
                i += 1;
                reference_path = Some(PathBuf::from(require_arg(args, i, "--reference")));
            }
            "--outdir" => {
                i += 1;
                outdir = PathBuf::from(require_arg(args, i, "--outdir"));
            }
            "--no-cyclic" => no_cyclic = true,
            "--no-reverse" => no_reverse = true,
            "--max-iter" => {
                i += 1;
                let v = require_arg(args, i, "--max-iter");
                max_iter = v.parse().map_err(|_| {
                    KnotError::DataParse(format!("bad --max-iter value '{v}' (expected integer)"))
                })?;
            }
            "--tol" => {
                i += 1;
                let v = require_arg(args, i, "--tol");
                tol = v.parse().map_err(|_| {
                    KnotError::DataParse(format!("bad --tol value '{v}' (expected float)"))
                })?;
            }
            other if !other.starts_with('-') => {
                return Err(KnotError::DataParse(format!(
                    "unexpected positional argument: {other}"
                )));
            }
            other => {
                return Err(KnotError::DataParse(format!("unknown flag '{other}'")));
            }
        }
        i += 1;
    }

    let selection_path =
        selection_path.ok_or_else(|| KnotError::DataParse("--selection is required".into()))?;
    let rows = read_knot_index(&selection_path)?;
    let row_by_frame = index_by_frame(&rows);

    let input_format = input_format.unwrap_or_else(|| detect_input_format(input_path));
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    let allow_reverse = !no_reverse;
    let cyclic_default = mode == Mode::FullChain;
    let cyclic = cyclic_default && !no_cyclic;

    match (input_format, mode) {
        (InputFormat::Xyz, Mode::FullChain) => run_full_chain(
            XyzFrameIter::new(reader),
            &row_by_frame,
            &target_type,
            arm_state,
            reference_path.as_deref(),
            &outdir,
            cyclic,
            allow_reverse,
            max_iter,
            tol,
        ),
        (InputFormat::Lammps, Mode::FullChain) => run_full_chain(
            LammpsFrameIter::new(reader),
            &row_by_frame,
            &target_type,
            arm_state,
            reference_path.as_deref(),
            &outdir,
            cyclic,
            allow_reverse,
            max_iter,
            tol,
        ),
        (InputFormat::Xyz, Mode::CoreByLength) => run_core_by_length(
            XyzFrameIter::new(reader),
            &row_by_frame,
            &target_type,
            arm_state,
            reference_path.as_deref(),
            &outdir,
            allow_reverse,
            max_iter,
            tol,
        ),
        (InputFormat::Lammps, Mode::CoreByLength) => run_core_by_length(
            LammpsFrameIter::new(reader),
            &row_by_frame,
            &target_type,
            arm_state,
            reference_path.as_deref(),
            &outdir,
            allow_reverse,
            max_iter,
            tol,
        ),
    }
}

fn run_full_chain<I>(
    iter: I,
    row_by_frame: &HashMap<usize, SelectionRow>,
    target_type: &str,
    arm_state: ArmState,
    reference_path: Option<&Path>,
    outdir: &Path,
    cyclic: bool,
    allow_reverse: bool,
    max_iter: usize,
    tol: f64,
) -> Result<()>
where
    I: Iterator<Item = Result<Vec<Point3>>>,
{
    let mut frames_two: Vec<Vec<Point3>> = Vec::new();
    let mut frames_three: Vec<Vec<Point3>> = Vec::new();
    let mut sel_two: Vec<usize> = Vec::new();
    let mut sel_three: Vec<usize> = Vec::new();

    for (idx, fr) in iter.enumerate() {
        let points = fr?;
        let row = match row_by_frame.get(&idx) {
            Some(r) => r,
            None => continue,
        };
        if row.knot_type != target_type {
            continue;
        }
        if !matches_arm(row, arm_state)? {
            continue;
        }

        match arm_state {
            ArmState::TwoArm => {
                frames_two.push(points);
                sel_two.push(idx);
            }
            ArmState::ThreeArm => {
                frames_three.push(points);
                sel_three.push(idx);
            }
            ArmState::Both => {
                let arm = row.arm_type.as_deref().ok_or_else(|| {
                    KnotError::DataParse("selection file missing arm_type column".into())
                })?;
                if arm == "two-arm" {
                    frames_two.push(points);
                    sel_two.push(idx);
                } else if arm == "three-arm" {
                    frames_three.push(points);
                    sel_three.push(idx);
                }
            }
        }
    }

    if arm_state == ArmState::TwoArm && frames_two.is_empty() {
        return Err(KnotError::DataParse("no frames selected (two-arm)".into()));
    }
    if arm_state == ArmState::ThreeArm && frames_three.is_empty() {
        return Err(KnotError::DataParse("no frames selected (three-arm)".into()));
    }
    if arm_state == ArmState::Both && frames_two.is_empty() && frames_three.is_empty() {
        return Err(KnotError::DataParse("no frames selected (two-arm or three-arm)".into()));
    }

    let reference = match reference_path {
        Some(p) => Some(load_reference(p)?),
        None => None,
    };

    match arm_state {
        ArmState::TwoArm => {
            let r = reference.as_deref().unwrap_or(&frames_two[0]);
            let (avg, results) = compute_average(&frames_two, r, cyclic, allow_reverse, max_iter, tol)?;
            write_average_outputs(outdir, &avg, &results, &sel_two)?;
        }
        ArmState::ThreeArm => {
            let r = reference.as_deref().unwrap_or(&frames_three[0]);
            let (avg, results) =
                compute_average(&frames_three, r, cyclic, allow_reverse, max_iter, tol)?;
            write_average_outputs(outdir, &avg, &results, &sel_three)?;
        }
        ArmState::Both => {
            if !frames_two.is_empty() {
                let sub = outdir.join("two-arm");
                let r = reference.as_deref().unwrap_or(&frames_two[0]);
                let (avg, results) =
                    compute_average(&frames_two, r, cyclic, allow_reverse, max_iter, tol)?;
                write_average_outputs(&sub, &avg, &results, &sel_two)?;
            }
            if !frames_three.is_empty() {
                let sub = outdir.join("three-arm");
                let r = reference.as_deref().unwrap_or(&frames_three[0]);
                let (avg, results) =
                    compute_average(&frames_three, r, cyclic, allow_reverse, max_iter, tol)?;
                write_average_outputs(&sub, &avg, &results, &sel_three)?;
            }
        }
    }
    Ok(())
}

#[derive(Default)]
struct GroupData {
    frames: Vec<Vec<Point3>>,
    selected: Vec<usize>,
}

fn run_core_by_length<I>(
    iter: I,
    row_by_frame: &HashMap<usize, SelectionRow>,
    target_type: &str,
    arm_state: ArmState,
    reference_path: Option<&Path>,
    outdir: &Path,
    allow_reverse: bool,
    max_iter: usize,
    tol: f64,
) -> Result<()>
where
    I: Iterator<Item = Result<Vec<Point3>>>,
{
    let reference = match reference_path {
        Some(p) => Some(load_reference(p)?),
        None => None,
    };

    let mut groups: HashMap<usize, GroupData> = HashMap::new();

    for (idx, fr) in iter.enumerate() {
        let points = fr?;
        let row = match row_by_frame.get(&idx) {
            Some(r) => r,
            None => continue,
        };
        if row.knot_type != target_type {
            continue;
        }
        if !matches_arm(row, arm_state)? {
            continue;
        }
        if row.knot_start < 0 || row.knot_end < 0 {
            continue;
        }
        let core = extract_core(&points, row)?;
        let len = core.len();
        if len == 0 {
            continue;
        }
        let g = groups.entry(len).or_default();
        g.frames.push(core);
        g.selected.push(idx);
    }

    if groups.is_empty() {
        return Err(KnotError::DataParse("no frames selected for core-by-length".into()));
    }

    std::fs::create_dir_all(outdir)?;

    let mut lens: Vec<usize> = groups.keys().copied().collect();
    lens.sort_unstable();

    for len in lens {
        let g = groups.get(&len).unwrap();
        if g.frames.is_empty() {
            continue;
        }
        let sub = outdir.join(format!("group_len_{:04}", len));
        let r = match reference.as_deref() {
            Some(rr) if rr.len() == len => rr,
            Some(_) => &g.frames[0],
            None => &g.frames[0],
        };
        let (avg, results) = compute_average(&g.frames, r, false, allow_reverse, max_iter, tol)?;
        write_average_outputs(&sub, &avg, &results, &g.selected)?;
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if let Err(e) = run_cli(&args[1..], &args[0]) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

