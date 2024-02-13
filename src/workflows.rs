
use std::collections::HashMap;

use crate::io::SampleOut;
use crate::recombination::RecombinationMapGenome;
use crate::variants::VCFCollection;
use log::info;
use rand::Rng;

// use rayon::prelude::*;

pub fn wrk_generate_offspring(sample: &SampleOut, genome_recomb_map: &RecombinationMapGenome, popvars: &VCFCollection, _denovo: &str, verbose: bool, contig_size: &HashMap<String, u64>, _n_threads: usize) -> () {
    
    if verbose {
        info!("Generating offspring for: {}", sample.name);
    }

    let cx_parent1 =
        genome_recomb_map.generate_genome_cx("parent1".to_string());
    let cx_parent2 = 
        genome_recomb_map.generate_genome_cx("parent2".to_string());

    let chr_vector = popvars.vcfs.keys().collect::<Vec<&String>>();
    for chr in chr_vector {

        let contig_size = contig_size.get(chr).unwrap().clone();
        // check if samples are in the vcf
        let vcf_obj = popvars.vcfs.get(chr).unwrap();

        if ! vcf_obj.samples.contains(&sample.parent1) {
            panic!("Sample {} not found in VCF", sample.parent1);
        }

        if !vcf_obj.samples.contains(&sample.parent2) {
            panic!("Sample {} not found in VCF", sample.parent2);
        }

        // these should be sorted
        let vcf_obj = popvars.vcfs.get(chr).unwrap();
        if verbose {
            info!("Chromosome: {}", chr);
        }
        let cx_chr_inst_p1 = cx_parent1.get(chr).unwrap();
        let cx_chr_inst_p2 = cx_parent2.get(chr).unwrap();

        let mut all_cx = Vec::new();
        all_cx.extend(cx_chr_inst_p1);
        all_cx.extend(cx_chr_inst_p2);

        // this is per chromosome
        // need to 
        // 1. choose a random haplotype for parent1
        let initial_haplotype_parent1: usize = rand::thread_rng().gen_range(0..2);
        // 2. choose a random haplotype for parent2
        let initial_haplotype_parent2: usize = rand::thread_rng().gen_range(0..2);

        // 3. combine cx from both parents, need to get a (hap1, hap2, position)
        let mut all_hap: Vec<(usize, usize, u64, u64)> = Vec::new();
        // note that this is 0-based
        let mut last_position: u64 = 0;

        let mut current_hap1 = initial_haplotype_parent1.clone();
        
        let mut current_hap2 = initial_haplotype_parent2.clone();

        if verbose {
            info!("Initial haplotype parent1: {}", current_hap1);
            info!("Initial haplotype parent2: {}", current_hap2);
        }

        // i am pretty sure if the length of cx is 0, this won't run, that is ok
        all_cx.iter().for_each(|(parent, crossover)| {
            let inst_position = crossover.position;
            all_hap.push((current_hap1.clone(),
                          current_hap2.clone(),
                          last_position,
                          inst_position));

            if parent == "parent1" {  // we should create a type for this?
                current_hap1 = if current_hap1 == 0 {1} else {0};
            } else {
                current_hap2 = if current_hap2 == 0 {1} else {0};
            }
        
            last_position = inst_position.clone();
        });

        // add last position, that is the end of the chromosome
        all_hap.push(
            (current_hap1.clone(),
             current_hap2.clone(),
             last_position,
             contig_size.clone())
        );

        // 4. iteratively print the positions where any of the parents has a variant
        // in that position 
        if verbose {
            info!("Haplotypes: {:?}", all_hap);
        }

        all_hap.iter().for_each(|(hap1, hap2, pos_from, pos_to)| {
                if verbose {
                    info!("Getting records for {}:{}-{}", chr, pos_from, pos_to);
                }
                vcf_obj.get_records_two_parents_from_to(
                    &sample.parent1,
                    &sample.parent2,
                    *hap1, 
                    *hap2,
                    chr.clone(),
                    *pos_from,
                    *pos_to,
                    verbose,
                    &sample.targetvcfout);
                });
    }
    ()
}



