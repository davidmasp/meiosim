
use crate::utils::ensure_directory_exists;

pub struct FamilyOut {
    pub prefix: String,
    pub samples: Vec<SampleOut>
}

pub struct SampleOut {
    pub name: String,
    pub parent1: String,
    pub parent2: String,
    pub targetvcfout: String,
    pub targetbedout: String,
}

impl SampleOut {
    fn new(prefix: &String, p1name: &String, p2name: &String, name: &String) -> Self {
        let targetvcfout = format!("{}/{}_{}_{}_meiosimvariants.txt", prefix, name, p1name, p2name);
        let targetbedout = format!("{}/{}_{}_{}_haplotypes.bed", prefix, name, p1name, p2name);
        Self {
            name: name.clone(),
            parent1: p1name.clone(),
            parent2: p2name.clone(),
            targetvcfout,
            targetbedout
        }
    }
}

impl FamilyOut {
    pub fn new(prefix: &String, parent1: &String, parent2: &String, number_of_sibs: u8, verbose: bool) -> Self {
        let mut samples = Vec::new();
        // check if the prefix exists as a folder name and if not create it
        ensure_directory_exists(&prefix, verbose).expect("Directory could not be created");
        for i in 0..number_of_sibs {
            let name = format!("sib{}", i);
            let sample = SampleOut::new(prefix,
                        parent1,
                        parent2,
                        &name);
            samples.push(sample);
        }
        Self {
            prefix: prefix.clone(),
            samples,
        }
    }
}

