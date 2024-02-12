

mod variants;
mod recombination;
mod utils;
use recombination::RecombinationMapGenome;

fn main() {

    // PARAMS
    let vcf_path = "debug/vcfcollections/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz";
    let path_with_rm = "debug/recombmaps";

    // VARIANTS
    let mut test = variants::VCF::new(vcf_path.to_string(), "chr22".to_string());
    let _ = test.read_samples();
    test.get_records( "HG02685".to_string(), 0).unwrap();

    // RECOM
    let genome_recomb_map = RecombinationMapGenome::from_path(path_with_rm, "map");
    let cx = genome_recomb_map.generate_genome_cx();

    for i in 0..cx.len() {
        println!("{}:{}", cx[i].seqname, cx[i].position);
    }

}
