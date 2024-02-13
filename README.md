
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
mkdir -p debug/vcfcollectionssmall2
files=$(ls debug/vcfcollectionssmall/*.vcf.gz)
for i in $files;
do
    ibase=$(basename $i)
    bcftools view -s NA21123,NA20752 $i | bcftools norm -m +snps | bcftools -o debug/vcfcollectionssmall2/$ibase
    bcftools index -t debug/vcfcollectionssmall2/$ibase
done
```


```
time ./target/release/meiosim \
    -r debug/recombmaps2/ \
    -v debug/vcfcollectionssmall2/   \
    -d ntt \
    -p NA21123 \
    -P NA20752 \
    --prefix testout \
    --seed 3 -f 5 --genome debug/hg38.genome --dwgsim
```

With the whole vcf the time is: 1m30s, with only the two parents, the time is: 10s

We also need to normalize the SNPs otherwise it's giving us problems
down the road. This essentially means that we are removing multiallelic
variants.

```
mutationfile=testout/sib0_NA21123_NA20752_meiosimvariants.txt
dwgsim -m ${mutationfile} -o 1 -c 0 -C 10 -R 0.0 smallgenome.fa outsib0
```

</details>
