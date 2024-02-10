

mod recombination;

fn main() {

    let rm_filename = "/home/dmas/src/DWGSIM/genetic_map_GRCh37_chr4.txt";
    let recomb_map_result = recombination::parse_to_recombination_map(rm_filename, true);

    let recomb_map = match recomb_map_result {
        Ok(rm) => {
            rm
        },
        Err(e) => {
            panic!("Error: {:?}", e);
        }
    };

    for i in 0..recomb_map.segments.len() {
        let next_segment_id = i + 1;
        if next_segment_id >= (recomb_map.segments.len()) {
            continue;
        }
        let next_segment = &recomb_map.segments[next_segment_id];
        let ncx = recomb_map.segments[i].sample_from_segment(next_segment);
        if ncx == 0 {
            continue;
        } else {
            let pos = recomb_map.segments[i].get_cx_position(next_segment, ncx);
            println!("CX in {:?}", pos);
        }
    }

}
