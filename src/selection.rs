use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

use crate::error::{KnotError, Result};

#[derive(Clone, Debug)]
pub struct SelectionRow {
    pub frame: usize,
    pub knot_type: String,
    pub knot_start: i32,
    pub knot_end: i32,
    pub knot_size: i32,
    pub arm_type: Option<String>,
}

fn parse_i32(field: &str, name: &str) -> Result<i32> {
    i32::from_str(field).map_err(|e| KnotError::DataParse(format!("bad {name} '{field}': {e}")))
}

fn parse_usize(field: &str, name: &str) -> Result<usize> {
    usize::from_str(field)
        .map_err(|e| KnotError::DataParse(format!("bad {name} '{field}': {e}")))
}

pub fn read_knot_index(path: &Path) -> Result<Vec<SelectionRow>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 5 {
            return Err(KnotError::DataParse(format!(
                "line {}: expected at least 5 columns, got {}",
                line_no + 1,
                parts.len()
            )));
        }
        let frame = parse_usize(parts[0], "frame")?;
        let knot_type = parts[1].to_string();
        let knot_start = parse_i32(parts[2], "knot_start")?;
        let knot_end = parse_i32(parts[3], "knot_end")?;
        let knot_size = parse_i32(parts[4], "knot_size")?;
        let arm_type = if parts.len() >= 6 {
            Some(parts[5].to_string())
        } else {
            None
        };
        rows.push(SelectionRow {
            frame,
            knot_type,
            knot_start,
            knot_end,
            knot_size,
            arm_type,
        });
    }
    if rows.is_empty() {
        return Err(KnotError::DataParse(format!(
            "no data rows found in selection file: {}",
            path.display()
        )));
    }
    Ok(rows)
}

pub fn index_by_frame(rows: &[SelectionRow]) -> HashMap<usize, SelectionRow> {
    let mut map = HashMap::with_capacity(rows.len());
    for r in rows {
        map.insert(r.frame, r.clone());
    }
    map
}

