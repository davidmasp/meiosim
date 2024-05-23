

use crate::utils::list_files_in_directory;
use crate::utils::capture_chromosome_from_file_name;
use crate::utils::from_vu8_to_string;

use core::panic;
use std::fs::File;
use std::io::Write;

use std::collections::HashMap;
use rust_htslib::bcf;
use bcf::Read;
use bcf::record::GenotypeAllele;

use log::warn;
use log::info;

use rand::Rng;
use rand::rngs::StdRng;


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
            output_writter: &mut File) -> (){
        // s.as_bytes()
        let fpath = self.file_path.clone();
        let mut bcf = bcf::IndexedReader::from_path(fpath)
            .expect("Cannot open the file");
        let header = bcf.header();
        let parent1_idx = header.sample_id(parent1.as_bytes()).expect("Parent1 not found");
        let parent2_idx = header.sample_id(parent2.as_bytes()).expect("Parent2 not found");

        let chr_id = header.name2rid(chromosome.as_bytes()).expect("Chromosome not found");
        bcf.fetch(chr_id, pos_from, Some(pos_to)).expect("Cannot fetch the region");

        let mut current_record = bcf.empty_record();

        while bcf.read(&mut current_record).is_some() {
            let pos = current_record.pos(); //  0-based position
            let pos1based = pos + 1;
            /*
                    I understand that for DWGSIM, the positions are 1-based
                    https://github.com/nh13/DWGSIM/blob/main/docs/03_Simulating_Reads.md#output-mutations-file
             */
            let alleles = from_vu8_to_string(current_record.alleles());

            let issnp = from_alleles_to_issnp(&alleles);
            if !issnp {
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
                    let line_out = compose_dwgsim_format(&chromosome, pos1based, alleles, parent1_gt_hapl, parent2_gt_hapl);
                    write!(output_writter, "{}", line_out).expect("Unable to write to file");
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

fn compose_dwgsim_format(chromosome: &String, pos1based: i64, mut alleles: Vec<String>, parent1_gt_hapl: &i32, parent2_gt_hapl: &i32) -> String {
    // see here https://github.com/nh13/DWGSIM/blob/main/docs/03_Simulating_Reads.md#output-mutations-file
    // I am unsure what the strand means here,
    // I am assuming all SNPs are strand 1 means from parent1 and 2 from parent 2.
    // the value 3 means that the SNP is homozygous, so that makes sense so far
    let mut string_out = String::new();
    if parent1_gt_hapl == &0 || parent2_gt_hapl == &0 && parent1_gt_hapl != parent2_gt_hapl{
        // this would be the case where the variant is heterozygous, as not 0|0 should be reported.
        // then we have to define the strand as 1 or 2.
        let strand: String;
        if parent1_gt_hapl == &0 {
            strand = "1".to_string();
        } else {
            strand = "2".to_string();
        }
        let ref_all = alleles[0].clone();
        // this is a bit stupid!
        let _ = alleles.sort();
        let iupac_code = get_iupac_representation(&alleles[0], &alleles[1]);
        string_out = format!("{}\t{}\t{}\t{}\t{}\n", chromosome, pos1based, ref_all, iupac_code, strand);
    } else if parent1_gt_hapl == &1 && parent2_gt_hapl == &1 {
        // this would be the case where the SNP is homozygous, strand does not matter     
        let ref_all = &alleles[0];
        let alt_all = &alleles[1];
        string_out = format!("{}\t{}\t{}\t{}\t3\n", chromosome, pos1based, ref_all, alt_all);
    }
    string_out
}

pub fn get_iupac_representation(x: &str, y: &str) -> String {
    // see here https://genome.ucsc.edu/goldenPath/help/iupac.html
    if x == "A" && y == "G" {
        "R".to_string()
    } else if x == "C" && y == "T" {
        "Y".to_string()
    } else if x == "C" && y == "G" {
        "S".to_string()
    } else if x == "A" && y == "T" {
        "W".to_string()
    } else if x == "G" && y == "T" {
        "K".to_string()
    } else if x == "A" && y == "C" {
        "M".to_string()
    } else {
        panic!("IUPAC code not found for {} and {}", x, y);
    }
}

pub fn flush_dnm_to_file(output_writter: &mut File, truepos_writter: &mut File, dnm: &String, verbose: bool, rng_dnm: &mut StdRng) -> () {
    // is this dumb?
    let dnm_reader1 = bcf::Reader::from_path(dnm)
            .expect("Cannot open the file");
    let header = dnm_reader1.header();

    let mut dnm_reader = bcf::Reader::from_path(dnm)
            .expect("Cannot open the file");

    dnm_reader.records().into_iter()
        .for_each(|x|{
            let record = x.unwrap();
            let pos = record.pos(); //  0-based position
            let pos1based = pos + 1;
            let chrom_u8 = Vec::from(header.rid2name(record.rid().unwrap()).unwrap());
            let chrom = String::from_utf8(chrom_u8).unwrap();
            let alleles = from_vu8_to_string(record.alleles());

            let issnp = from_alleles_to_issnp(&alleles);
            if !issnp {
                if verbose {
                    warn!("Skipping record at {}:{} because it's not SNP/SNV", chrom, pos);
                }
                return;
            }

            // I am forcing them to be "heterozygous" but
            // the parent of choice is "random"
            
            let gts = if rng_dnm.gen_range(0..2) == 0 { (0, 1) } else { (1, 0) };
        
            let line_out = compose_dwgsim_format(
                        &chrom,
                        pos1based,
                        alleles,
                        &gts.0, 
                        &gts.1);
            write!(output_writter, "{}", line_out).expect("Unable to write to file");
            write!(truepos_writter, "{}\t{}\n", chrom, pos1based).expect("Unable to write to file");
        });
        ()
}

pub fn flush_vcf_to_file(output_writter: &mut File, vcfname: &String, verbose: bool) -> () {
    // is this dumb?
    let vcf_reader1 = bcf::Reader::from_path(vcfname)
            .expect("Cannot open the file");
    let header = vcf_reader1.header();

    let mut vcf_reader = bcf::Reader::from_path(vcfname)
            .expect("Cannot open the file");

    vcf_reader.records().into_iter()
        .for_each(|x|{
            let record = x.unwrap();
            let pos = record.pos(); //  0-based position
            let pos1based = pos + 1;
            let chrom_u8 = Vec::from(header.rid2name(record.rid().unwrap()).unwrap());
            let chrom = String::from_utf8(chrom_u8).unwrap();
            let alleles = from_vu8_to_string(record.alleles());

            let issnp = from_alleles_to_issnp(&alleles);
            if !issnp {
                if verbose {
                    warn!("Skipping record at {}:{} because it's not SNP/SNV", chrom, pos);
                }
                return;
            }
            // currently we are assuming the vcf is uni-sample here as I
            // have it implemented it like this in the pipeline,
            // this can be improved by adding a specific sample.
            // I think however this way is easier to parall.
            let gts_raw: bcf::record::Genotype = record.genotypes().unwrap().get(0);
            let gt1_ga: &GenotypeAllele = gts_raw.get(0).unwrap();
            let gt2_ga : &GenotypeAllele = gts_raw.get(1).unwrap();

            // the first record should always be unphased
            // https://docs.rs/rust-htslib/latest/rust_htslib/bcf/record/struct.Genotypes.html#implementations
            let gt1_int = match gt1_ga {
                GenotypeAllele::Unphased(value) => value,
                GenotypeAllele::Phased(_value) => panic!("problem with phasing??"),
                GenotypeAllele::UnphasedMissing => panic!("Missing GT error"),
                GenotypeAllele::PhasedMissing => panic!("Phased missing GT error"),
            };

            let gt2_int = match gt2_ga {
                GenotypeAllele::Unphased(_value) => panic!("problem with phasing??"),
                GenotypeAllele::Phased(value) => value,
                GenotypeAllele::UnphasedMissing => panic!("Missing GT error"),
                GenotypeAllele::PhasedMissing => panic!("Phased missing GT error"),
            };

            let line_out = compose_dwgsim_format(
                        &chrom,
                        pos1based,
                        alleles,
                        &gt1_int, 
                        &gt2_int);
            write!(output_writter, "{}", line_out).expect("Unable to write to file");
        });
        ()
}

fn from_alleles_to_issnp(alleles: &Vec<String>) -> bool {
    let alleles_len: Vec<usize>= alleles.iter()
        .map(|allele| {allele.len()})
        .collect();
    let any_greater_than_one = alleles_len.iter().any(|&len| len > 1);

    if alleles.len() != 2{
        false
    } else if any_greater_than_one{
        false
    } else {
        true
    }
}

