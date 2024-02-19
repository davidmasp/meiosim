
use crate::utils::list_files_in_directory;

use rand::Rng;
use rand::rngs::StdRng;

use rand_distr::{Poisson, Distribution};

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::collections::HashMap;
use std::cmp;

#[derive(Clone, Debug, Eq, Default, Hash)]
pub struct Crossover {
    pub seqname: String,
    pub position: u64
}

impl PartialOrd for Crossover {
    fn partial_cmp(&self, other: &Crossover) -> Option<cmp::Ordering> {
        Some(other.cmp(self))
    }
}

impl Ord for Crossover {
    fn cmp(&self, other: &Crossover) -> cmp::Ordering {
        self.position.cmp(&other.position)
    }
}

impl PartialEq for Crossover {
    fn eq(&self, other: &Crossover) -> bool {
        self.position == other.position && self.seqname == other.seqname
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Default)]
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
    pub fn sample_from_segment(&self, next_segment: &RecombinationSegment, rng: &mut StdRng) -> u8 {
        if self.centimorgan == next_segment.centimorgan {
            return 0;
        }
        let cmdistance = next_segment.centimorgan - self.centimorgan;
        //println!("dist: {}", cmdistance);
        let cxlambda = cmdistance * 0.01;
        //println!("lambda: {}", cxlambda);
        let poi = Poisson::new(cxlambda).unwrap();
        let recombination_count = poi.sample(rng);
        let recombination_count_u8: u8 = recombination_count as u8;
        recombination_count_u8
    }
    pub fn get_cx_position(&self, next_segment: &RecombinationSegment, recombination_count: u8, rnd: &mut StdRng) -> Vec<u64> {
        let mut positions = Vec::new();
        for _ in 0..recombination_count {
            let pos_inst = rnd.gen_range(self.position..next_segment.position);
            positions.push(pos_inst);
        }
        // here we need to order the positions, also added a test to check
        positions.sort();
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
    pub fn generate_cx(&self, rng_cx: &mut StdRng) -> Vec<Crossover>{
        let mut vec_out = Vec::new();
        for i in 0..self.segments.len() {
            let next_segment_id = i + 1;
            if next_segment_id >= (self.segments.len()) {
                continue;
            }
            let next_segment = &self.segments[next_segment_id];
            let ncx = self.segments[i].sample_from_segment(next_segment, rng_cx);
            if ncx == 0 {
                continue;
            } else {
                let pos_cx = self.segments[i]
                    .get_cx_position(next_segment,
                                     ncx,
                                     rng_cx);
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
    pub fn from_path(path: &str, extension: &str, recom_header: bool) -> Self {
        let rm_filenames = list_files_in_directory(path, extension).unwrap();
        let recombination_maps = rm_filenames.iter()
            .map(|rm_filename|{
                let recomb_map_result = RecombinationMap::parse_to_recombination_map(rm_filename, recom_header);
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
    pub fn generate_genome_cx(&self, parentid: String, rng_cx: &mut StdRng) -> HashMap<String, Vec<(String, Crossover)>> {
        let mut map_out = HashMap::new();
        for recombination_map in &self.recombination_maps {
            let mut cx = recombination_map.generate_cx(rng_cx);
            cx.sort();
            let cxout: Vec<(String, Crossover)> = cx.iter().cloned().map(|c| (parentid.clone(), c.clone())).collect();
            map_out.insert(recombination_map.seqname.clone(), cxout);
        }
        map_out
    }
}


#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use super::RecombinationSegment;

    #[test]
    fn test_get_cx_position() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let seg1: RecombinationSegment = RecombinationSegment::new("chr1".to_string(), 0, 0.0);
        let seg2 = RecombinationSegment::new("chr1".to_string(), 100000, 10.0);
        let positions = seg1.get_cx_position(&seg2, 5, &mut rng);
        println!("{:?}", positions);

        // Check that the positions vector is not empty
        assert!(!positions.is_empty(), "Positions vector is empty");

        // Check that the positions are ordered
        let mut last_position = positions[0];
        for &position in positions.iter().skip(1) {
            assert!(position >= last_position, "Positions are not ordered");
            last_position = position;
        }
    }
}

