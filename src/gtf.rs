// Copyright (C) 2025 [Thomas Steimlé]
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Represents a gene with its name and genomic coordinates.
#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub start: u64,
    pub end: u64,
}

/// Loads a GTF file and returns a HashMap mapping chromosome names to vectors of Gene structs.
/// Only "gene" features are extracted from the GTF file.
///
/// # Arguments
/// * `path` - Path to the GTF file
///
/// # Returns
/// * A HashMap where keys are chromosome names and values are vectors of Gene structs
pub fn load_gtf(path: &Path) -> anyhow::Result<HashMap<String, Vec<Gene>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut genes_map: HashMap<String, Vec<Gene>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        // Skip comment lines and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        // Use iterator approach to check field count without collecting all fields
        let mut fields_iter = line.split('\t');
        let chrom = match fields_iter.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        // Skip source (field 1)
        fields_iter.next();
        let feature_type = match fields_iter.next() {
            Some(f) => f,
            None => continue,
        };

        // Only process "gene" features
        if feature_type != "gene" {
            continue;
        }

        let start_str = match fields_iter.next() {
            Some(s) => s,
            None => continue,
        };
        let end_str = match fields_iter.next() {
            Some(e) => e,
            None => continue,
        };

        // Parse start and end, skip invalid entries
        let start: u64 = match start_str.parse() {
            Ok(s) => s,
            Err(_) => continue, // Skip entries with invalid start coordinates
        };
        let end: u64 = match end_str.parse() {
            Ok(e) => e,
            Err(_) => continue, // Skip entries with invalid end coordinates
        };

        // Skip score (field 5), strand (field 6), frame (field 7)
        for _ in 0..3 {
            if fields_iter.next().is_none() {
                continue;
            }
        }

        let attributes = match fields_iter.next() {
            Some(a) => a,
            None => continue,
        };

        // Extract gene_name or gene_id from attributes
        let gene_name = extract_gene_name(attributes);

        if let Some(name) = gene_name {
            let gene = Gene { name, start, end };
            genes_map.entry(chrom).or_default().push(gene);
        }
    }

    Ok(genes_map)
}

/// Extracts gene_name or gene_id from GTF attributes column.
/// Prefers gene_name over gene_id if both are present.
fn extract_gene_name(attributes: &str) -> Option<String> {
    let mut gene_name: Option<String> = None;
    let mut gene_id: Option<String> = None;

    for attr in attributes.split(';') {
        let attr = attr.trim();
        if attr.is_empty() {
            continue;
        }

        // GTF format: key "value" or key value
        let parts: Vec<&str> = attr.splitn(2, ' ').collect();
        if parts.len() < 2 {
            continue;
        }

        let key = parts[0].trim();
        let value = parts[1].trim().trim_matches('"');

        match key {
            "gene_name" => gene_name = Some(value.to_string()),
            "gene_id" => gene_id = Some(value.to_string()),
            _ => {}
        }
    }

    // Prefer gene_name over gene_id
    gene_name.or(gene_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_gene_name_with_gene_name() {
        let attrs = r#"gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_type "transcribed_unprocessed_pseudogene";"#;
        assert_eq!(extract_gene_name(attrs), Some("DDX11L1".to_string()));
    }

    #[test]
    fn test_extract_gene_name_with_only_gene_id() {
        let attrs = r#"gene_id "ENSG00000223972"; gene_type "transcribed_unprocessed_pseudogene";"#;
        assert_eq!(extract_gene_name(attrs), Some("ENSG00000223972".to_string()));
    }

    #[test]
    fn test_extract_gene_name_empty() {
        let attrs = "";
        assert_eq!(extract_gene_name(attrs), None);
    }
}
