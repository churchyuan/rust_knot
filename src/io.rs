use std::io::{BufRead, Write};

use crate::error::{KnotError, Result};
use crate::point::Point3;

/// Lazy frame-by-frame iterator over an XYZ stream.
/// Reads one frame at a time without loading the entire file.
pub struct XyzFrameIter<R> {
    reader: R,
    line_buf: String,
    done: bool,
}

impl<R: BufRead> XyzFrameIter<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
            done: false,
        }
    }
}

impl<R: BufRead> Iterator for XyzFrameIter<R> {
    type Item = Result<Vec<Point3>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Find next non-empty header line (point count)
        let n: usize = loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => {
                    self.done = true;
                    return None;
                }
                Ok(_) => {}
                Err(e) => return Some(Err(KnotError::Io(e))),
            }
            let trimmed = self.line_buf.trim();
            if trimmed.is_empty() {
                continue;
            }
            match trimmed.parse() {
                Ok(v) => break v,
                Err(e) => {
                    return Some(Err(KnotError::DataParse(format!(
                        "invalid point count: {e}"
                    ))))
                }
            }
        };

        // Skip comment line
        self.line_buf.clear();
        match self.reader.read_line(&mut self.line_buf) {
            Ok(0) => {
                return Some(Err(KnotError::DataParse(
                    "unexpected EOF after XYZ point count".into(),
                )))
            }
            Err(e) => return Some(Err(KnotError::Io(e))),
            _ => {}
        }

        // Read n data lines
        let mut points = Vec::with_capacity(n);
        for i in 0..n {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => {
                    return Some(Err(KnotError::DataParse(format!(
                        "unexpected EOF in XYZ data line {}",
                        i + 1
                    ))))
                }
                Err(e) => return Some(Err(KnotError::Io(e))),
                _ => {}
            }
            let parts: Vec<&str> = self.line_buf.split_whitespace().collect();
            if parts.len() < 4 {
                return Some(Err(KnotError::DataParse(format!(
                    "line {}: expected 'type x y z', got '{}'",
                    i + 3,
                    self.line_buf.trim()
                ))));
            }
            let coords: std::result::Result<[f64; 3], _> = (|| {
                Ok([
                    parts[1]
                        .parse()
                        .map_err(|e| format!("bad x at line {}: {e}", i + 3))?,
                    parts[2]
                        .parse()
                        .map_err(|e| format!("bad y at line {}: {e}", i + 3))?,
                    parts[3]
                        .parse()
                        .map_err(|e| format!("bad z at line {}: {e}", i + 3))?,
                ])
            })();
            match coords {
                Ok(c) => points.push(c),
                Err(msg) => return Some(Err(KnotError::DataParse(msg))),
            }
        }

        Some(Ok(points))
    }
}

/// Lazy frame-by-frame iterator over a LAMMPS dump stream.
/// Expects the standard repeated block:
/// ITEM: TIMESTEP
/// ...
/// ITEM: NUMBER OF ATOMS
/// N
/// ITEM: BOX BOUNDS ...
/// ...
/// ITEM: ATOMS id type x y z
/// ...
pub struct LammpsFrameIter<R> {
    reader: R,
    line_buf: String,
    done: bool,
}

impl<R: BufRead> LammpsFrameIter<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
            done: false,
        }
    }
}

impl<R: BufRead> Iterator for LammpsFrameIter<R> {
    type Item = Result<Vec<Point3>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match parse_lammps_frame(&mut self.reader, &mut self.line_buf) {
            Ok(Some(points)) => Some(Ok(points)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(err) => {
                self.done = true;
                Some(Err(err))
            }
        }
    }
}

/// Read points from XYZ format.
/// Format:
/// ```text
/// N
/// (comment line)
/// atom_type  x  y  z
/// ...
/// ```
pub fn read_data_xyz<R: BufRead>(reader: &mut R) -> Result<Vec<Point3>> {
    let mut line = String::new();

    // Read N
    reader.read_line(&mut line).map_err(KnotError::Io)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| KnotError::DataParse(format!("invalid point count: {e}")))?;

    // Skip comment line
    line.clear();
    reader.read_line(&mut line).map_err(KnotError::Io)?;

    let mut points = Vec::with_capacity(n);
    for i in 0..n {
        line.clear();
        reader.read_line(&mut line).map_err(KnotError::Io)?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(KnotError::DataParse(format!(
                "line {}: expected 'type x y z', got '{}'",
                i + 3,
                line.trim()
            )));
        }
        let x: f64 = parts[1]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad x at line {}: {e}", i + 3)))?;
        let y: f64 = parts[2]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad y at line {}: {e}", i + 3)))?;
        let z: f64 = parts[3]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad z at line {}: {e}", i + 3)))?;
        points.push([x, y, z]);
    }

    Ok(points)
}

/// Read all frames from an XYZ stream into memory.
/// For large files prefer `XyzFrameIter` with batched processing.
pub fn read_data_xyz_frames<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<Point3>>> {
    // Delegate to iterator — avoids duplicated parsing logic
    let mut frames = Vec::new();
    // Take ownership via a wrapper that implements BufRead
    // Since we have &mut R, construct iterator over it
    let iter = XyzFrameIter::new(reader);
    for result in iter {
        frames.push(result?);
    }
    Ok(frames)
}

/// Read points from LAMMPS dump format.
/// Expects: 3 header lines, then N, then 5 more header lines, then data.
pub fn read_data_lammps<R: BufRead>(reader: &mut R) -> Result<Vec<Point3>> {
    let mut line = String::new();
    parse_lammps_frame(reader, &mut line)?
        .ok_or_else(|| KnotError::DataParse("no LAMMPS frame found".into()))
}

/// Read all frames from a LAMMPS dump stream into memory.
pub fn read_data_lammps_frames<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<Point3>>> {
    let mut frames = Vec::new();
    let iter = LammpsFrameIter::new(reader);
    for result in iter {
        frames.push(result?);
    }
    Ok(frames)
}

/// Write points in XYZ format.
pub fn write_data_xyz<W: Write>(points: &[Point3], writer: &mut W) -> Result<()> {
    writeln!(writer, "{}", points.len()).map_err(KnotError::Io)?;
    writeln!(writer).map_err(KnotError::Io)?;
    for p in points {
        writeln!(writer, "1\t{:10.5}\t{:10.5}\t{:10.5}", p[0], p[1], p[2])
            .map_err(KnotError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_xyz() {
        let data = "3\ncomment\n1 1.0 2.0 3.0\n1 4.0 5.0 6.0\n1 7.0 8.0 9.0\n";
        let mut reader = Cursor::new(data);
        let points = read_data_xyz(&mut reader).unwrap();
        assert_eq!(points.len(), 3);
        assert_eq!(points[0], [1.0, 2.0, 3.0]);
        assert_eq!(points[2], [7.0, 8.0, 9.0]);
    }

    #[test]
    fn test_write_read_roundtrip() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut buf = Vec::new();
        write_data_xyz(&points, &mut buf).unwrap();

        let mut reader = Cursor::new(buf);
        let read_back = read_data_xyz(&mut reader).unwrap();
        assert_eq!(read_back.len(), 2);
        // Check approximate equality due to formatting
        for (a, b) in points.iter().zip(read_back.iter()) {
            for d in 0..3 {
                assert!((a[d] - b[d]).abs() < 0.001);
            }
        }
    }

    #[test]
    fn test_read_xyz_frames() {
        let data = concat!(
            "2\n",
            "frame0\n",
            "1 1.0 2.0 3.0\n",
            "1 4.0 5.0 6.0\n",
            "3\n",
            "frame1\n",
            "1 7.0 8.0 9.0\n",
            "1 10.0 11.0 12.0\n",
            "1 13.0 14.0 15.0\n"
        );
        let mut reader = Cursor::new(data);
        let frames = read_data_xyz_frames(&mut reader).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].len(), 2);
        assert_eq!(frames[1].len(), 3);
        assert_eq!(frames[0][0], [1.0, 2.0, 3.0]);
        assert_eq!(frames[1][2], [13.0, 14.0, 15.0]);
    }

    #[test]
    fn test_read_xyz_frames_with_blank_lines_and_empty_input() {
        let data = "\n\n1\nframe0\n1 1.0 2.0 3.0\n\n";
        let mut reader = Cursor::new(data);
        let frames = read_data_xyz_frames(&mut reader).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0][0], [1.0, 2.0, 3.0]);

        let mut empty_reader = Cursor::new("");
        let empty_frames = read_data_xyz_frames(&mut empty_reader).unwrap();
        assert!(empty_frames.is_empty());
    }

    #[test]
    fn test_read_lammps_frames() {
        let data = concat!(
            "ITEM: TIMESTEP\n",
            "0\n",
            "ITEM: NUMBER OF ATOMS\n",
            "2\n",
            "ITEM: BOX BOUNDS ff ff ff\n",
            "0 1\n",
            "0 1\n",
            "0 1\n",
            "ITEM: ATOMS id type x y z\n",
            "1 1 1.0 2.0 3.0\n",
            "2 1 4.0 5.0 6.0\n",
            "ITEM: TIMESTEP\n",
            "1\n",
            "ITEM: NUMBER OF ATOMS\n",
            "1\n",
            "ITEM: BOX BOUNDS ff ff ff\n",
            "0 1\n",
            "0 1\n",
            "0 1\n",
            "ITEM: ATOMS id type x y z\n",
            "1 1 7.0 8.0 9.0\n"
        );
        let mut reader = Cursor::new(data);
        let frames = read_data_lammps_frames(&mut reader).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0][0], [1.0, 2.0, 3.0]);
        assert_eq!(frames[0][1], [4.0, 5.0, 6.0]);
        assert_eq!(frames[1][0], [7.0, 8.0, 9.0]);
    }
}

fn parse_lammps_frame<R: BufRead>(
    reader: &mut R,
    line: &mut String,
) -> Result<Option<Vec<Point3>>> {
    line.clear();
    loop {
        match reader.read_line(line).map_err(KnotError::Io)? {
            0 => return Ok(None),
            _ if !line.trim().is_empty() => break,
            _ => line.clear(),
        }
    }

    if line.trim() != "ITEM: TIMESTEP" {
        return Err(KnotError::DataParse(format!(
            "expected 'ITEM: TIMESTEP', got '{}'",
            line.trim()
        )));
    }

    line.clear();
    if reader.read_line(line).map_err(KnotError::Io)? == 0 {
        return Err(KnotError::DataParse(
            "unexpected EOF after LAMMPS timestep header".into(),
        ));
    }

    line.clear();
    if reader.read_line(line).map_err(KnotError::Io)? == 0 {
        return Err(KnotError::DataParse(
            "unexpected EOF before atom-count header".into(),
        ));
    }
    if line.trim() != "ITEM: NUMBER OF ATOMS" {
        return Err(KnotError::DataParse(format!(
            "expected 'ITEM: NUMBER OF ATOMS', got '{}'",
            line.trim()
        )));
    }

    line.clear();
    if reader.read_line(line).map_err(KnotError::Io)? == 0 {
        return Err(KnotError::DataParse("unexpected EOF at atom count".into()));
    }
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| KnotError::DataParse(format!("invalid atom count: {e}")))?;

    line.clear();
    if reader.read_line(line).map_err(KnotError::Io)? == 0 {
        return Err(KnotError::DataParse(
            "unexpected EOF before box-bounds header".into(),
        ));
    }
    if !line.trim().starts_with("ITEM: BOX BOUNDS") {
        return Err(KnotError::DataParse(format!(
            "expected 'ITEM: BOX BOUNDS...', got '{}'",
            line.trim()
        )));
    }

    for _ in 0..3 {
        line.clear();
        if reader.read_line(line).map_err(KnotError::Io)? == 0 {
            return Err(KnotError::DataParse("unexpected EOF in box bounds".into()));
        }
    }

    line.clear();
    if reader.read_line(line).map_err(KnotError::Io)? == 0 {
        return Err(KnotError::DataParse(
            "unexpected EOF before atom data header".into(),
        ));
    }
    if !line.trim().starts_with("ITEM: ATOMS") {
        return Err(KnotError::DataParse(format!(
            "expected 'ITEM: ATOMS...', got '{}'",
            line.trim()
        )));
    }

    let mut points = Vec::with_capacity(n);
    for _ in 0..n {
        line.clear();
        if reader.read_line(line).map_err(KnotError::Io)? == 0 {
            return Err(KnotError::DataParse(
                "unexpected EOF in atom coordinates".into(),
            ));
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            return Err(KnotError::DataParse(format!(
                "LAMMPS data line too short: '{}'",
                line.trim()
            )));
        }
        let x: f64 = parts[2]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad x: {e}")))?;
        let y: f64 = parts[3]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad y: {e}")))?;
        let z: f64 = parts[4]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad z: {e}")))?;
        points.push([x, y, z]);
    }

    Ok(Some(points))
}
