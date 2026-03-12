use std::collections::HashMap;
use std::io::BufRead;

use crate::error::{KnotError, Result};
use crate::polynomial::{parse_polynomial, Polynomial};

/// Table mapping Alexander polynomials to knot type names.
///
/// Since the Alexander polynomial is not a complete invariant, one polynomial
/// may correspond to multiple knot types (e.g. `8_3` and `10_1` share
/// `4-9*t+4*t^2`). The table stores all candidates per polynomial.
///
/// - `lookup()` returns the simplest candidate (lowest crossing number) — backward compatible.
/// - `lookup_all()` returns all candidates for advanced disambiguation.
#[derive(Debug)]
pub struct AlexanderTable {
    table: HashMap<Polynomial, Vec<String>>,
}

/// Parse "3_1" -> (3, 1) for sorting by crossing number.
fn knot_sort_key(name: &str) -> (u32, u32) {
    let parts: Vec<&str> = name.split('_').collect();
    if parts.len() == 2 {
        if let (Ok(c), Ok(i)) = (parts[0].parse(), parts[1].parse()) {
            return (c, i);
        }
    }
    // unknot "1" or unparseable — sort first
    (0, 0)
}

impl AlexanderTable {
    /// Parse the table from a reader (file or string).
    /// Format: each line is `knot_name\tpolynomial_string` or space-separated.
    ///
    /// In strict mode (default), malformed lines cause an error with line number.
    /// In non-strict mode, malformed lines are silently skipped.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        Self::from_reader_impl(reader, true)
    }

    /// Parse with explicit strict flag.
    pub fn from_reader_lenient<R: BufRead>(reader: R) -> Result<Self> {
        Self::from_reader_impl(reader, false)
    }

    fn from_reader_impl<R: BufRead>(reader: R, strict: bool) -> Result<Self> {
        let mut table: HashMap<Polynomial, Vec<String>> = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line.map_err(KnotError::Io)?;
            let line = line.trim().to_string();
            if line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.splitn(2, |c: char| c == '\t' || c == ' ').collect();
            if parts.len() < 2 {
                if strict {
                    return Err(KnotError::DataParse(format!(
                        "line {}: expected 'knot_name polynomial', got '{line}'",
                        line_num + 1
                    )));
                }
                continue;
            }

            let knot_name = parts[0].trim().to_string();
            let poly_str = parts[1].trim();
            let poly = parse_polynomial(poly_str)?;

            // Insert under both poly and -poly (Alexander polynomial defined up to sign)
            table
                .entry(poly.clone())
                .or_default()
                .push(knot_name.clone());

            let neg_poly = poly.negate();
            if neg_poly != poly {
                table.entry(neg_poly).or_default().push(knot_name);
            }
        }

        // Sort each candidate list by crossing number (simplest first)
        for candidates in table.values_mut() {
            candidates.sort_by_key(|n| knot_sort_key(n));
            candidates.dedup();
        }

        Ok(AlexanderTable { table })
    }

    /// Load from a file path.
    pub fn from_file(path: &str) -> Result<Self> {
        let file = std::fs::File::open(path).map_err(KnotError::Io)?;
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
    }

    /// Look up the simplest (lowest crossing number) knot type for a polynomial.
    /// Tries both the polynomial and its negation.
    pub fn lookup(&self, poly: &Polynomial) -> Option<&str> {
        if let Some(names) = self.table.get(poly) {
            if let Some(first) = names.first() {
                return Some(first.as_str());
            }
        }
        let neg = poly.negate();
        self.table
            .get(&neg)
            .and_then(|names| names.first())
            .map(|s| s.as_str())
    }

    /// Look up all candidate knot types for a polynomial, sorted by crossing number.
    /// Returns an empty slice if not found.
    pub fn lookup_all(&self, poly: &Polynomial) -> &[String] {
        if let Some(names) = self.table.get(poly) {
            if !names.is_empty() {
                return names;
            }
        }
        let neg = poly.negate();
        self.table
            .get(&neg)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Number of distinct polynomials in the table (including negations).
    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_table() {
        let data = "1\t1\n3_1\t-1+t-t^2\n4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        assert!(!table.is_empty());

        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));

        let neg_trefoil = trefoil.negate();
        assert_eq!(table.lookup(&neg_trefoil), Some("3_1"));

        let unknot = parse_polynomial("1").unwrap();
        assert_eq!(table.lookup(&unknot), Some("1"));
    }

    #[test]
    fn test_figure_eight() {
        let data = "4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(table.lookup(&poly), Some("4_1"));
    }

    #[test]
    fn test_ambiguous_polynomial() {
        // Two knots share the same polynomial
        let data = "8_3\t4-9*t+4*t^2\n10_1\t4-9*t+4*t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("4-9*t+4*t^2").unwrap();

        // lookup returns simplest (lowest crossing number)
        assert_eq!(table.lookup(&poly), Some("8_3"));

        // lookup_all returns both, sorted by crossing number
        let all = table.lookup_all(&poly);
        assert_eq!(all, &["8_3", "10_1"]);
    }

    #[test]
    fn test_lookup_all_unique() {
        let data = "3_1\t-1+t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("-1+t-t^2").unwrap();

        assert_eq!(table.lookup_all(&poly), &["3_1"]);
    }

    #[test]
    fn test_lookup_all_not_found() {
        let data = "3_1\t-1+t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("999+t").unwrap();

        assert!(table.lookup_all(&poly).is_empty());
    }

    #[test]
    fn test_strict_mode_rejects_bad_line() {
        let data = "3_1\t-1+t-t^2\nbadline\n4_1\t-1+3*t-t^2\n";
        let result = AlexanderTable::from_reader(data.as_bytes());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("line 2"));
    }

    #[test]
    fn test_lenient_mode_skips_bad_line() {
        let data = "3_1\t-1+t-t^2\nbadline\n4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader_lenient(data.as_bytes()).unwrap();
        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));
        let fig8 = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(table.lookup(&fig8), Some("4_1"));
    }
}
