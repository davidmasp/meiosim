
use crate::io::SampleOut;
use crate::recombination::RecombinationMapGenome;
use crate::variants::VCFCollection;
use crate::variants as vars;
use rand::Rng;

pub fn wrk_generate_offspring(sample: &SampleOut, genome_recomb_map: &RecombinationMapGenome, popvars: &VCFCollection, denovo: &str, verbose: bool) -> () {
    if verbose {
        println!("Generating offspring for: {}", sample.name);
    }

    let cx_parent1 = genome_recomb_map.generate_genome_cx();
    let cx_parent2 = genome_recomb_map.generate_genome_cx();

    let chr_vector = popvars.vcfs.keys().collect::<Vec<&String>>();
    for chr in chr_vector {
        // check if samples are in the vcf
        let vcf_obj = popvars.vcfs.get(chr).unwrap();

        if ! vcf_obj.samples.contains(&sample.parent1) {
            panic!("Sample {} not found in VCF", sample.parent1);
        }

        if !vcf_obj.samples.contains(&sample.parent2) {
            panic!("Sample {} not found in VCF", sample.parent2);
        }

        // this is per chromosome
        // need to 
        // 1. choose a random haplotype for parent1
        let haplotype_parent1: usize = rand::thread_rng().gen_range(0..2);
        // 2. choose a random haplotype for parent2
        let haplotype_parent2: usize = rand::thread_rng().gen_range(0..2);
        // 3. iteratively print the positions where any of the parents has a variant in that position 
        let vcf_obj = popvars.vcfs.get(chr).unwrap();
        vcf_obj.get_records_two_parents(
                &sample.parent1,
                &sample.parent2,
                 haplotype_parent1, 
                 haplotype_parent2);     
    }
    ()
}



