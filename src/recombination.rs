

use rand::Rng;
use rand_distr::{Poisson, Distribution};

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub struct RecombinationSegment {
    pub seqname: String,
    pub position: u64,
    pub rate: f64,
    pub centimorgan: f64,
}

pub struct RecombinationMap {
    pub seqname: String,
    pub segments: Vec<RecombinationSegment>,
}

pub struct RecombinationMapCollection {
    pub recombination_maps: Vec<RecombinationMap>,
}

impl RecombinationSegment {
    fn new(seqname: String, position: u64, rate: f64, centimorgan: f64) -> Self {
        Self {
            seqname,
            position,
            rate,
            centimorgan,
        }
    }
    pub fn create_initial_segment(seqname: String) -> Self {
        Self {
            seqname,
            position: 0,
            rate: 0.0,
            centimorgan: 0.0,
        }
    }
    pub fn sample_from_segment(&self, next_segment: &RecombinationSegment) -> u8 {
        // is this okay? 
        // this could come from outside the function but I think there
        // is implications for parall. processing

        if self.centimorgan == next_segment.centimorgan {
            return 0;
        }

        let mut rng = rand::thread_rng();

        let cmdistance = next_segment.centimorgan - self.centimorgan;
        //println!("dist: {}", cmdistance);
        let cxlambda = cmdistance * 0.01;
        //println!("lambda: {}", cxlambda);
        let poi = Poisson::new(cxlambda).unwrap();
        let recombination_count = poi.sample(&mut rng);
        let recombination_count_u8: u8 = recombination_count as u8;
        recombination_count_u8
    }
    pub fn get_cx_position(&self, next_segment: &RecombinationSegment, recombination_count: u8) -> Vec<u64> {
        let mut rng = rand::thread_rng();
        let mut positions = Vec::new();

        for _ in 0..recombination_count {
            let pos_inst = rng.gen_range(self.position..next_segment.position);
            positions.push(pos_inst);
        }
        positions
    }
}

impl RecombinationMap {
    fn new(seqname: String, segments: Vec<RecombinationSegment>) -> Self {
        Self {
            seqname,
            segments,
        }
    }
}

pub fn parse_to_recombination_map(file_path: &str, has_header: bool) -> io::Result<RecombinationMap> {
    
    // note that one file always needs to be one chromosome!!!
    let input_file = File::open(&Path::new(file_path))?;
    let reader = io::BufReader::new(input_file);
    let mut lines = reader.lines();

    if has_header {
        let _ = lines.next(); // Skip header line
    }

    let mut current_seqname = String::new();
    let mut segments = Vec::new();

    for line in lines {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 4 {
                let seqname = parts[0].to_string();
                if current_seqname != seqname && !segments.is_empty() {
                    panic!("Multiple chromosomes in one file are not supported");
                }
                let position = parts[1].parse::<u64>().unwrap();
                let rate = parts[2].parse::<f64>().unwrap();
                let centimorgan = parts[3].parse::<f64>().unwrap();
                current_seqname = seqname.clone();
                segments.push(RecombinationSegment::new(seqname, position, rate, centimorgan));
            }
        }
    }

    let recombination_map = RecombinationMap::new(current_seqname, segments);

    Ok(recombination_map)
}

