
// Tests for the gtf_parser library
#[cfg(test)]
mod tests {
    use gtf_gene_structure::GTF::Gtf;
    use gtf_gene_structure::ExonIterator::ExonIterator;
    use gtf_gene_structure::GTF::QueryErrors;
    use gtf_gene_structure::Gene::RegionStatus;

    #[test]
    fn test_parse_gtf() {
        let mut gtf = Gtf::new();
        let result = gtf.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path
        
        // Check if the parsing was successful
        assert!(result.is_ok());

        // Add your assertions here to verify the contents of the Gtf
        // For example:
        assert_eq!(gtf.chromosomes.len(), 2);

        // Verify the names of the chromosomes (replace "chr1" and "chr2" with your expected names)
        let expected_chromosomes = vec!["Y".to_string(), "Y_mapped_Scaffold_9_D1573".to_string() ]; // Replace with actual expected names
        let mut actual_chromosomes: Vec<String> = gtf.chromosomes.keys().cloned().collect();
        actual_chromosomes.sort();
        assert_eq!(actual_chromosomes, expected_chromosomes);
    }

    #[test]
    fn test_slice_gtf() {
        let mut gtf = Gtf::new(); // Use lowercase `gtf` for the variable
        let _ = gtf.parse_gtf("tests/data/test.gtf"); // Use a valid test GTF file path

        // Assuming we have added some test data here for the gtf

        // Perform a slice operation
        let sliced = match gtf.slice_Gtf("Y", 3265434, 3276434) {
            Ok(slice) => slice,
            Err(e) => panic!("{:?}", e),
        };
        
        // Check the expected result
        assert_eq!(sliced.len(), 4);
    }

    #[test]
    fn test_query() {
        let mut gtf = Gtf::new();

        // Add exons on chromosome "Y"
        gtf.add_exon("gene1_id", "gene1_name", 20, 40, "Y".to_string(), true);
        gtf.add_exon("gene2_id", "gene2_name", 50, 80, "Y".to_string(), true);
        gtf.add_exon("gene3_id", "gene3_name", 100, 200, "Y".to_string(), true);


        // Perform a search on chromosome "Y" with specific start-end positions
        let mut interator = ExonIterator::new("test");
        gtf.init_search( "Y", 10, &mut interator );

        let mut result_y = match gtf.query("Y", 10, 15, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::BeforeGene );

        result_y = match gtf.query("Y", 15, 25, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );

        result_y = match gtf.query("Y", 20, 30, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::InsideExon );

        result_y = match gtf.query("Y", 10, 50, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );

        result_y = match gtf.query("Y", 10, 50, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::SpanningBoundary );
        assert_eq!(result_y.0, ("gene1_id".to_string(), "gene1_name".to_string() ));

        result_y = match gtf.query("Y", 60, 80, &mut interator ) {
            Ok(slice) => slice,
            Err(e) => panic!("This should not throw an error!"),
        };
        assert_eq!(result_y.1, RegionStatus::InsideExon );
        assert_eq!(result_y.0, ("gene2_id".to_string(), "gene2_name".to_string() ));

        let _ = match gtf.query("Y", 6000, 80000, &mut interator ){
            Ok(slice) => panic!("This should fail here") ,
            Err(e) => {
                assert_eq!(e, QueryErrors::OutOfGenes )
            },
        };
    }

}
