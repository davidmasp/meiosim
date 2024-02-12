
# Meiosim

Simulate meiosis from a panel of population vcfs (aka 1kG).

## Input

### Recombination maps

Recombination maps per chromosome need to be inputed. You can get
recombination maps for Human (GRCh38) in the
[beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip) docs.

The maps need to be in
[PLINK map format](https://zzz.bwh.harvard.edu/plink/data.shtml#map).

### VCF collections from 1KG

The vcf collections (aka 1 multisample-VCF per chromosome)
are available from [1000genomes EBI ftp site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).


<details>

```
./target/debug/meiosim -r sdfjsd -v dsjhfkds -d dsjhfkjs -p SAMPLE1 -P SAMPLE2 --prefix test --seed 2 -f 10
```

</details>
