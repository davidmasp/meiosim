

use crate::utils::list_files_in_directory;
use crate::utils::capture_chromosome_from_file_name;

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::collections::HashMap;


use clap::Error;
use noodles::vcf::record::Chromosome;
use noodles::vcf::record::Position;
// noodles imports
// use indexmap::IndexSet;
use noodles::vcf as vcf;
use indexmap::IndexSet;
use noodles::vcf::Record;
use noodles::vcf::record::genotypes::keys::key::GENOTYPE as GenoKey;

pub struct VCF {
    pub file_path: String,
    pub seqname: String,
    pub samples: IndexSet<String>,
}

impl VCF {
    pub fn new(file_path: String, seqname: String) -> VCF {
        VCF {
            file_path,
            seqname,
            samples: IndexSet::new(),
        }
    }
    pub fn get_reader(&self) -> vcf::Reader<Box<dyn BufRead>> {
        let reader = vcf::reader::Builder::default().build_from_path(&self.file_path).expect("msg");
        reader
    }
    pub fn read_samples(&mut self) -> io::Result<()> {
        let mut reader = self.get_reader();
        let header = reader.read_header()?;
        let samples = header.sample_names();
        //println!("Samples: {:?}", samples);
        self.samples = samples.clone();
        Ok(())
    }
    fn get_record(record: &Record, haplotype: &usize, sample_idx: &usize) -> Option<Record> {
        // this is a point of optimization
        let sample_idx_clone = sample_idx.clone();
        let haplotype_clone = haplotype.clone();
        let genotype_of_interest = record.genotypes()
            .get_index(sample_idx_clone)
            .expect("Sample not found");
        let short_geno_option = genotype_of_interest
            .get(&GenoKey)
            .expect("Genotype not found");
        let short_geno_value: String = match short_geno_option {
            Some(short_geno) => {
                short_geno.clone().to_string()
            },
            None => {
                panic!("Genotypes: None");
            }
        };
        let fields: Vec<&str> = short_geno_value.split("|").collect();
        // mmmmh
        if fields.len() != 2 {
            panic!("Genotypes: {:?} (?not phased)", short_geno_value);
        }
        let genotype_of_selected_haplotype = fields[haplotype_clone].to_string();
        if genotype_of_selected_haplotype == "0" {
            None
        } else {
            Some(record.clone())
        }
    }
    pub fn _get_records(&self, sample_of_interest: String, haplotype: usize) -> io::Result<Vec<Record>> {
        let mut reader = self.get_reader();
        let header = reader.read_header()?;
        // here i need to get the index of my sample of interest
        let sample_idx = self.samples.get_index_of(&sample_of_interest).expect("dfsdfs");
        let mut records_out = Vec::new();
        for result in reader.records(&header) {
            let record: Record = result.expect("dhfkhsf");

            let record_option = VCF::get_record(&record, &haplotype, &sample_idx);

            match record_option {
                Some(record) => {
                    records_out.push(record);
                },
                None => {
                    continue;
                }
            }
        }
        Ok(records_out)
    }
    pub fn get_records_two_parents(&self, parent1: &String, parent2: &String, haplotype1: usize, haplotype2: usize) -> () {
        let mut reader = self.get_reader();
        let header = reader.read_header().expect("problem reading header");

        let parent1_idx = self.samples.get_index_of(parent1).expect("dfsdfs");
        let parent2_idx = self.samples.get_index_of(parent2).expect("dfsdfs");

        for result in reader.records(&header) {
            let record: Record = result.expect("dhfkhsf");

            let record_parent1_option = VCF::get_record(&record, &haplotype1, &parent1_idx);
            let record_parent2_option = VCF::get_record(&record, &haplotype2, &parent2_idx);

            if let (None, None) = (&record_parent1_option, &record_parent2_option) {
                continue;
            } else if let (genop1, None) = (&record_parent1_option, &record_parent2_option) {
                println!("{}\t{}\t{}|{}", record.chromosome(), record.position(), "1", "0");
            } else if let (None, genop2) = (&record_parent1_option, &record_parent2_option) {
                println!("{}\t{}\t{}|{}", record.chromosome(), record.position(), "0", "1");
            } else {
                println!("{}\t{}\t{}|{}", record.chromosome(), record.position(), "1", "1");
            }
        }
        ()
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
            let mut vcf: VCF = VCF::new(file, seqname);
            let _ = vcf.read_samples();
            vcfs.insert(seqname2, vcf);
        }

        if vcfs.len() == 0 {
            panic!("No VCF files found in the folder");
        }

        VCFCollection {
            vcfs,
        }
    }

    pub fn get_vcf(&self, seqname: &str) -> Option<&VCF> {
        self.vcfs.get(seqname)
    }
}
