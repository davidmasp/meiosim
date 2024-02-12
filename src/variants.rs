


use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::collections::HashMap;

// noodles imports
// use indexmap::IndexSet;
use noodles::vcf as vcf;
use indexmap::IndexSet;
use noodles::vcf::Record;
use noodles::vcf::record::genotypes::keys::key::GENOTYPE as GenoKey;

pub struct VCF {
    file_path: String,
    seqname: String,
    samples: IndexSet<String>,
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
    pub fn get_records(&self, sample_of_interest: String, haplotype: usize) -> io::Result<Vec<Record>> {
        let mut reader = self.get_reader();
        let header = reader.read_header()?;
        // here i need to get the index of my sample of interest
        let sample_idx = self.samples.get_index_of(&sample_of_interest).expect("dfsdfs");
        let mut records_out = Vec::new();
        for result in reader.records(&header) {
            let record = result.expect("dhfkhsf");
            let genotype_of_interest = record.genotypes().get_index(sample_idx).expect("fksdjh");
            let short_geno_option = genotype_of_interest.get(&GenoKey).expect("kdhfkjsh");
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
                panic!("Genotypes: {:?}", short_geno_value);
            }
            let genotype_of_selected_haplotype = fields[haplotype].to_string();
            if genotype_of_selected_haplotype == "0" {
                continue;
            } else {
                println!("Genotypes: {:?}", short_geno_value);
                records_out.push(record);
            }
        }
        Ok(records_out)
    }
}

pub struct VCFCollection {
    vcfs: HashMap<String, VCF>,
}

impl VCFCollection {
    pub fn new() -> VCFCollection {
        VCFCollection {
            vcfs: HashMap::new(),
        }
    }

    pub fn add_vcf(&mut self, vcf: VCF) {
        self.vcfs.insert(vcf.seqname.clone(), vcf);
    }

    pub fn get_vcf(&self, seqname: &str) -> Option<&VCF> {
        self.vcfs.get(seqname)
    }
}

