
use crate::utils::list_files_in_directory;

use rand::Rng;
use rand_distr::{Poisson, Distribution};

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::collections::HashMap;

pub struct Crossover {
    pub seqname: String,
    pub position: u64
}

pub struct RecombinationSegment {
    pub seqname: String,
    pub position: u64,
    pub centimorgan: f64,
}

pub struct RecombinationMap {
    pub seqname: String,
    pub segments: Vec<RecombinationSegment>,
}

pub struct RecombinationMapGenome {
    pub recombination_maps: Vec<RecombinationMap>,
}

impl RecombinationSegment {
    fn new(seqname: String, position: u64, centimorgan: f64) -> Self {
        Self {
            seqname,
            position,
            centimorgan,
        }
    }
    pub fn _create_initial_segment(seqname: String) -> Self {
        Self {
            seqname,
            position: 0,
            centimorgan: 0.0,
        }
    }
    pub fn sample_from_segment(&self, next_segment: &RecombinationSegment) -> u8 {
        // is this okay? 
        // this could come from outside the function but I think there
        // is implications for parall. processing
        let mut rng = rand::thread_rng();

        if self.centimorgan == next_segment.centimorgan {
            return 0;
        }
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
    pub fn parse_to_recombination_map(file_path: &str, has_header: bool) -> io::Result<Self> {
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
                let parts: Vec<&str> = line.split(' ').collect();
                if parts.len() == 4 {
                    let seqname = parts[0].to_string();
                    if current_seqname != seqname && !segments.is_empty() {
                        panic!("Multiple chromosomes in one file are not supported");
                    }
                    let position = parts[3].parse::<u64>().unwrap();
                    let centimorgan = parts[2].parse::<f64>().unwrap();
                    current_seqname = seqname.clone();
                    segments.push(RecombinationSegment::new(seqname, position, centimorgan));
                } else {
                    panic!("Invalid number of columns in recombination map file");
                }
            }
        }

        let recombination_map = RecombinationMap::new(current_seqname, segments);

        Ok(recombination_map)
    }
    pub fn generate_cx(&self) -> Vec<Crossover>{
        let mut vec_out = Vec::new();
        for i in 0..self.segments.len() {
            let next_segment_id = i + 1;
            if next_segment_id >= (self.segments.len()) {
                continue;
            }
            let next_segment = &self.segments[next_segment_id];
            let ncx = self.segments[i].sample_from_segment(next_segment);
            if ncx == 0 {
                continue;
            } else {
                let pos_cx = self.segments[i].get_cx_position(next_segment, ncx);
                let cx_inst = Crossover {
                    seqname: self.segments[i].seqname.clone(),
                    position: pos_cx[0],
                };
                vec_out.push(cx_inst);
            }
        }
        vec_out
    }
}

impl RecombinationMapGenome {
    pub fn from_path(path: &str, extension: &str) -> Self {
        let rm_filenames = list_files_in_directory(path, extension).unwrap();
        let recombination_maps = rm_filenames.iter()
            .map(|rm_filename|{
                let recomb_map_result = RecombinationMap::parse_to_recombination_map(rm_filename, true);
                let recomb_map = match recomb_map_result {
                    Ok(rm) => {
                        rm
                    },
                    Err(e) => {
                        panic!("Error: {:?}", e);
                    }
                };
                recomb_map
            })
            .collect();
        Self {
            recombination_maps,
        }
    }
    pub fn generate_genome_cx(&self) -> HashMap<String, Vec<Crossover>> {
        let mut map_out = HashMap::new();
        for recombination_map in &self.recombination_maps {
            let cx = recombination_map.generate_cx();
            map_out.insert(recombination_map.seqname.clone(), cx);
        }
        map_out
    }
}
