

use crate::utils::list_files_in_directory;
use crate::utils::capture_chromosome_from_file_name;
use crate::utils::from_vu8_to_string;

use std::fs::File;
use std::io::Write;

use std::collections::HashMap;
use rust_htslib::bcf;
use bcf::Read;
use bcf::record::GenotypeAllele;

use log::warn;
use log::info;

pub struct VCF {
    pub file_path: String,
    pub seqname: String,
    pub samples: Vec<String>,
}

impl VCF {
    pub fn new(file_path: String, seqname: String) -> VCF {
        let file_path2 = file_path.clone();
        let bcf = bcf::Reader::from_path(file_path).expect("Error opening VCF file");
        let header = bcf.header();
        let sample_names_vu8: Vec<& [u8]> = header.samples();
        let sample_names = from_vu8_to_string(sample_names_vu8);

        VCF {
            file_path: file_path2,
            seqname,
            samples: sample_names,
        }
    }

    pub fn get_records_two_parents_from_to(&self,
            parent1: &String, 
            parent2: &String, 
            hap1: usize, 
            hap2: usize, 
            chromosome: String, 
            pos_from: u64, 
            pos_to: u64, 
            verbose: bool,
            output_file: &String) -> (){
        // s.as_bytes()
        let fpath = self.file_path.clone();
        let mut bcf = bcf::IndexedReader::from_path(fpath)
            .expect("Cannot open the file");
        let header = bcf.header();
        let parent1_idx = header.sample_id(parent1.as_bytes()).expect("Parent1 not found");
        let parent2_idx = header.sample_id(parent2.as_bytes()).expect("Parent2 not found");

        let chr_id = header.name2rid(chromosome.as_bytes()).expect("Chromosome not found");
        bcf.fetch(chr_id, pos_from, Some(pos_to)).expect("Cannot fetch the region");

        // open the file
        let mut outputfile = File::create(output_file).expect("Unable to create file");

        let mut current_record = bcf.empty_record();

        while bcf.read(&mut current_record).is_some() {
            let pos = current_record.pos(); //  0-based position
            let pos1based = pos + 1;
            /*
                    I understand that for DWGSIM, the positions are 1-based
                    https://github.com/nh13/DWGSIM/blob/main/docs/03_Simulating_Reads.md#output-mutations-file
             */
            let alleles = from_vu8_to_string(current_record.alleles());
            let alleles_len: Vec<usize>= alleles.iter()
                .map(|allele| {allele.len()})
                .collect();
            let any_greater_than_one = alleles_len.iter().any(|&len| len > 1);
            if alleles.len() != 2{
                if verbose {
                    warn!("Skipping record at {}:{} because it's multiallelic", chromosome, pos);
                }
                continue; 
            } else if any_greater_than_one{
                if verbose {
                    warn!("Skipping record at {}:{} because it's not SNP/SNV", chromosome, pos);
                }
                continue;
            }

            let genotypes = current_record.genotypes().expect("Error reading genotypes");

            let parent1_gt = genotypes.get(parent1_idx);
            let parent2_gt = genotypes.get(parent2_idx);

            let parent1_gt_hapl = extract_value(parent1_gt.get(hap1).unwrap());
            let parent2_gt_hapl = extract_value(parent2_gt.get(hap2).unwrap());

            if parent1_gt_hapl.is_none() || parent2_gt_hapl.is_none() {
                if verbose {
                    warn!("Skipping record at {}:{} because one of the parents has missing genotype", chromosome, pos);
                }
                continue;
            } else {
                let parent1_gt_hapl = parent1_gt_hapl.unwrap();
                let parent2_gt_hapl = parent2_gt_hapl.unwrap();
                if parent1_gt_hapl == &(0 as i32) && parent2_gt_hapl == &(0 as i32) {
                    continue;
                } else {
                    write!(outputfile, "{}\t{}\t{:?}\t{}|{}\n", chromosome, pos1based, alleles, parent1_gt_hapl, parent2_gt_hapl).expect("Unable to write to file");
                }
            }
        }

    }
}

pub struct VCFCollection {
    pub vcfs: HashMap<String, VCF>,
}

impl VCFCollection {
    pub fn from_path(path: &str, extension: &str, verbose: bool) -> VCFCollection {
        let list_of_files = list_files_in_directory(path, extension).unwrap();

        let mut vcfs = HashMap::new();
        for file in list_of_files {
            let seqname = capture_chromosome_from_file_name(&file, verbose);
            let seqname2 = seqname.clone();
            let vcf: VCF = VCF::new(file, seqname);
            vcfs.insert(seqname2, vcf);
        }

        if vcfs.len() == 0 {
            panic!("No VCF files found in the folder");
        } else {
            if verbose {
                info!("Total VCFs found: {:?}", vcfs.len());
            }
        }

        VCFCollection {
            vcfs,
        }
    }
}


fn extract_value(allele: &GenotypeAllele) -> Option<&i32> {
    match allele {
        GenotypeAllele::Unphased(value) => Some(value),
        GenotypeAllele::Phased(value) => Some(value),
        GenotypeAllele::UnphasedMissing => None,
        GenotypeAllele::PhasedMissing => None,
    }
}
