
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
        Some(self.cmp(other))
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
            let ncx_usize = ncx as usize;
            if ncx == 0 {
                continue;
            } else {
                let pos_cx = self.segments[i]
                    .get_cx_position(next_segment,
                                     ncx,
                                     rng_cx);
                for posid in 0..ncx_usize {
                    let cx_inst = Crossover {
                        seqname: self.segments[i].seqname.clone(),
                        position: pos_cx[posid],
                    };
                    vec_out.push(cx_inst);
                }
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
            // these should not be needed actually
            cx.sort();
            let cxout: Vec<(String, Crossover)> = cx.iter().cloned().map(|c| (parentid.clone(), c.clone())).collect();
            map_out.insert(recombination_map.seqname.clone(), cxout);
        }
        map_out
    }
}

pub fn generate_simple_recombination(csize: &HashMap<String, u64>, parentid: String, rng_cx: &mut StdRng, cxcount: u8) -> HashMap<String, Vec<(String, Crossover)>> {
    let mut map_out = HashMap::new();
    for (seqname, size) in csize.iter() {
        // size here means the size of the chromosome
        // and seqname should be the names of the chromosome
        // we can generate two recombination segments, the cM values do not
        // matter as they won't be used.
        // we cannot use 0 here because then i think it could generate a 0 (although unlikely)
        let segment1 = RecombinationSegment::new(seqname.clone(), 1, 0.0);
        let segment2 = RecombinationSegment::new(seqname.clone(), *size, 0.0);
        let positions = segment1.get_cx_position(&segment2, cxcount, rng_cx);
        
        let cxout: Vec<(String, Crossover)> = positions.iter().map(|pos| {
            (parentid.clone(), Crossover {
                seqname: seqname.clone(),
                position: *pos,
            })
        }).collect();
        map_out.insert(seqname.clone(), cxout);
    }
    map_out
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rand::Rng;
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

    #[test]
    fn test_crossover_order() {
        let mut rng: StdRng = StdRng::seed_from_u64(44);
        let seg1 = RecombinationSegment::new("chr1".to_string(), 0, 0.0);
        let seg2 = RecombinationSegment::new("chr1".to_string(), 100000, 500.0);
        let seg3 = RecombinationSegment::new("chr1".to_string(), 200000, 1000.0);
        // tested with this seed this test sould give more than 1 crossover

        let recom_map = super::RecombinationMap::new("chr1".to_string(), vec![seg1, seg2, seg3]);
        let crossovers = recom_map.generate_cx(&mut rng);
        // Check that the crossovers vector is not empty
        assert!(!crossovers.is_empty(), "Crossovers vector is empty");
        assert!(!crossovers.len() > 1, "Crossovers vector is equal to 1");
        println!("{:?}", crossovers);
        // Check that the crossovers are in descending order
        let mut last_position = crossovers[0].position;
        for crossover in crossovers.iter().skip(1) {
            assert!(crossover.position > last_position, "Crossovers are not in ascending order");
            last_position = crossover.position;
        }
    }

    #[test]
    fn test_rng_0() {
        let mut rng: StdRng = StdRng::seed_from_u64(44);
        let pos1 = &rng.gen_range(0..1);
        let pos1prime = &rng.gen_range(0..1);
        let pos2 = &rng.gen_range(1..2);
        let pos3 = &rng.gen_range(2..3);
        assert!(*pos1 == 0, "pos1 is not 1");
        assert!(*pos2 == 1, "pos1 is not 2");
        assert!(*pos3 == 2, "pos1 is not 3");
        assert!(*pos1prime == 0, "pos1prime is not 0");
    }
}

