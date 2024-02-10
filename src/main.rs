

mod recombination;
mod utils;
use recombination::RecombinationMapGenome;

fn main() {

    let path_with_rm = "debug/recombmaps";
    let genome_recomb_map = RecombinationMapGenome::from_path(path_with_rm, "map");
    let cx = genome_recomb_map.generate_genome_cx();

    for i in 0..cx.len() {
        println!("{}:{}", cx[i].seqname, cx[i].position);
    }

}
