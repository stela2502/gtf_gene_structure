use std::collections::HashMap;
use std::cmp::Ordering;
use std::error::Error;

use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use std::fmt;


#[derive(Debug, Clone, Copy)]
struct Exon {
    start: u32,
    end: u32,
}


#[derive(Debug, Clone)]
struct Gene {
    pub gene_id: String,
    start: u32,
    end: u32,
    sens_orientation: bool,
    exons: Vec<Exon>, // New field to store exons
}


#[derive(Debug, PartialEq)]
pub enum RegionStatus {
    InsideExon,
    InsideIntron,
    SpanningBoundary,
    OutsideGene,
}


#[derive(Debug)]
struct ExonIterator<'a> {
    exons: &'a [Exon],
    current_index: usize,
}

impl<'a> ExonIterator<'a> {
    // Constructor to create a new ExonIterator with a reference to the exons
    fn new(exons: &'a [Exon]) -> Self {
        ExonIterator {
            exons,
            current_index: 0,
        }
    }

    // Method to check if the given entry is in an exon
    fn is_in_exon(&mut self, position: u32) -> Option<&'a Exon> {
        while self.current_index < self.exons.len() {
            let exon = &self.exons[self.current_index];

            // If the position is within the exon, return it
            if position >= exon.start && position <= exon.end {
                return Some(exon);
            }

            // Move to the next exon
            self.current_index += 1;
        }

        // If no exon is found, return None
        None
    }

    // Reset the iterator if needed
    fn reset(&mut self) {
        self.current_index = 0;
    }
}

impl Gene {
    // Constructor for Gene
    pub fn new(gene_id: String, start: u32, end: u32, sens_orientation: bool) -> Self {
        Gene {
            gene_id,
            start,
            end,
            sens_orientation,
            exons: Vec::new(), // Initialize the exons vector
        }
    }

    // Method to add an exon to the gene
    pub fn add_exon(&mut self, start: u32, end: u32) {
        let exon = Exon { start, end };
        self.exons.push(exon); // Add the exon to the gene's exon list

        // Optionally, update the gene's start and end to reflect the new exon
        if start < self.start {
            self.start = start;
        }
        if end > self.end {
            self.end = end;
        }
    }

    // Method to check the region
    pub fn check_region(&self, region_start: u32, region_end: u32) -> RegionStatus {
        let mut is_within_exon = false;
        let mut is_within_intron = false;
        let mut spans_boundary = false;

        for exon in &self.exons {
            // Check if the region is entirely within an exon
            if region_start >= exon.start && region_end <= exon.end {
                is_within_exon = true;
            }

            // Check if the region is entirely outside the exons (this indicates it might be in an intron)
            if (region_start < exon.start && region_end < exon.start) || (region_start > exon.end && region_end > exon.end) {
                is_within_intron = true;
            }

            // Check if the region spans the exon-intron boundary
            if (region_start < exon.start && region_end > exon.start) || (region_start < exon.end && region_end > exon.end) {
                spans_boundary = true;
            }
        }

        if is_within_exon {
            RegionStatus::InsideExon
        } else if is_within_intron {
            RegionStatus::InsideIntron
        } else if spans_boundary {
            RegionStatus::SpanningBoundary
        } else {
            RegionStatus::OutsideGene
        }
    }
}


#[derive(Debug)]
struct Genome {
    chromosomes: HashMap<String, Vec<Gene>>, // Group genes and exons by chromosome
}

// Implement the Display trait for Genome
impl fmt::Display for Genome {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let num_chromosomes = self.chromosomes.len();
        writeln!(f, "Number of chromosomes: {}", num_chromosomes)?;

        for (chromosome, gene_exons) in &self.chromosomes {
            let num_genes = gene_exons.len();
            writeln!(f, "Chromosome {}: {} genes", chromosome, num_genes)?;
        }

        Ok(())
    }
}


impl Genome {
    fn new() -> Self {
        Genome {
            chromosomes: HashMap::new(),
        }
    }

    pub fn add_exon(&mut self, gene_id: String, start: u32, end: u32, chromosome: String, sens_orientation: bool) {
        // Get the vector of genes for the specified chromosome
        let chromosome_genes = self.chromosomes.entry(chromosome.clone()).or_insert(Vec::new());

        // Check if the gene already exists in the vector
        if let Some(gene) = chromosome_genes.iter_mut().find(|g| g.gene_id == gene_id) {
            // If it exists, add the exon to the existing gene
            gene.add_exon(start, end);
        } else {
            // If the gene does not exist, create a new gene and add the exon
            let mut new_gene = Gene::new(gene_id, start, end, sens_orientation);
            new_gene.add_exon(start, end);
            chromosome_genes.push(new_gene);
        }
    }

    fn sort_genes(&mut self) {
        for gene_exons in self.chromosomes.values_mut() {
            // Sort by gene start position, then exon start position
            gene_exons.sort_by(|a, b | {
                let gene_cmp = a.start.cmp(&b.start);
                if gene_cmp == Ordering::Equal {
                    a.gene_id.cmp(&b.gene_id) // Sort by gene ID to maintain consistency
                } else {
                    gene_cmp
                }
            });
        }
    }


    pub fn slice_genome(&self, chromosome: &str, start: u32, end: u32) -> Result<Vec<Gene>, &str> {
        // Check if the chromosome exists in the genome
        if let Some(genes) = self.chromosomes.get(chromosome) {
            // Filter genes based on their start and end positions
            Ok( genes.iter()
                .filter(|gene| gene.end > start && gene.start < end) // Ensure overlapping genes
                .cloned() // Clone each gene for the resulting vector
                .collect() // Collect results into a vector
            )
        } else {
            // If the chromosome does not exist, return an empty vector
            println!("Sorry chromosome {} is not present here", chromosome);
            Err("Sorry this chromosome is not present here")
        }
    }

    // Function to parse the GTF file and populate the Genome structure
    pub fn parse_gtf(&mut self, file_path: &str ) -> Result<(), Box<dyn Error>> {
        let path = Path::new(file_path);
        let file = File::open(&path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue; // Skip comment lines
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue; // Skip malformed lines
            }

            let chromosome = fields[0].to_string();
            let feature_type = fields[2];
            let orientation = if fields[7] == "+" { true } else { false };

            if feature_type == "exon" {
                let start: u32 = fields[3].parse()?;
                let end: u32 = fields[4].parse()?;

                // Extract gene ID from the last column attributes
                let attributes = fields[8];
                let gene_id = attributes
                    .split(';')
                    .find_map(|attr| {
                        let mut parts = attr.trim().split_whitespace();
                        if let Some(key) = parts.next() {
                            if key == "gene_id" {
                                parts.next() // Get the value after gene_id
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    })
                    .unwrap_or("unknown")
                    .to_string();

                self.add_exon(gene_id, start, end, chromosome, orientation);
            }
        }

        println!("I have read this:\n{}", self);
        Ok(())
    }
    
}



// Tests for the gtf_parser library
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gtf() {
        let mut genome = Genome::new();
        let result = genome.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path
        
        // Check if the parsing was successful
        assert!(result.is_ok());

        // Add your assertions here to verify the contents of the genome
        // For example:
        assert_eq!(genome.chromosomes.len(), 2);

        // Verify the names of the chromosomes (replace "chr1" and "chr2" with your expected names)
        let expected_chromosomes = vec!["Y".to_string(), "Y_mapped_Scaffold_9_D1573".to_string() ]; // Replace with actual expected names
        let mut actual_chromosomes: Vec<String> = genome.chromosomes.keys().cloned().collect();
        actual_chromosomes.sort();
        assert_eq!(actual_chromosomes, expected_chromosomes);
    }

    #[test]
    fn test_slice_genome() {
        let mut genome = Genome::new();
        let _ = genome.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path
        // Assuming we have added some test data here for the genome

        // Perform a slice operation
        let sliced = match genome.slice_genome("Y", 3265434, 3276434){
            Ok(slice) => slice,
            Err(e) => panic!("{:?}",e),
        };
        
        // Check the expected result
        assert_eq!(sliced.len(), 4);
    }
}

