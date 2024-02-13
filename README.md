
# Meiosim

Simulate meiosis from a panel of population vcfs (aka 1kG).

## Input

### Recombination maps

Recombination maps per chromosome need to be inputed. You can get
recombination maps for Human (GRCh38) in the
[beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip) docs.

The maps need to be in
[PLINK map format](https://zzz.bwh.harvard.edu/plink/data.shtml#map).

In the debug case, we also need to change the chromosome names like so:

```
files=$(ls debug/recombmaps/*.map)
mkdir -p debug/recombmaps2
for i in $files;
do
    ibase=$(basename $i)
    awk '{print("chr"$0)}' $i > debug/recombmaps2/${ibase}
done
```

### VCF collections from 1KG

The vcf collections (aka 1 multisample-VCF per chromosome)
are available from [1000genomes EBI ftp site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).


<details>

```
time ./target/debug/meiosim \
    -r debug/recombmaps2/ \
    -v debug/vcfcollectionssmall/   \
    -d ntt \
    -p NA21123 \
    -P NA20752 \
    --prefix testout \
    --seed 3 -f 10 --genome debug/hg38.genome
```

</details>
