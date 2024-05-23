
# Meiosim

Simulate meiosis from a panel of population vcfs (aka 1kG).

## Input

### Recombination maps

Recombination maps per chromosome need to be inputed. You can get
recombination maps for Human (GRCh38) in the
[beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip) docs.

The maps need to be in
[PLINK map format](https://zzz.bwh.harvard.edu/plink/data.shtml#map).

### VCF population data collections from 1KG

The vcf collections (aka 1 multisample-VCF per chromosome)
are available from [1000genomes EBI ftp site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).

### *de novo* variant collections

The DECODE dataset contain +1k trios with an average DNM of 67 mutations/trio.
They can be used (almost) directly with meiosim to introduce DNM by selecting
one sample randomly, you can download this data from their [paper](https://www.nature.com/articles/nature24018#Sec28).

<details>

## Prepare testing data

### Population variants

Download, chr21 and chr22 from [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).
This will take a while.

Once this is done, you can downsample the vcfs to ~10Mb with:

```	
vcfname="1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
mkdir -p debug/vcfcollectionssmall
bcftools view -r chr21:0-10000000 -o debug/vcfcollectionssmall/chr21.vcf.gz ${vcfname}
```

Then to speed up the program you can also select the samples that you are going to
use, it's also important to normalize the snps other-wise dwgsim will complain.

```
mkdir -p debug/vcfcollectionssmall2
files=$(ls debug/vcfcollectionssmall/*.vcf.gz)
for i in $files;
do
    ibase=$(basename $i)
    bcftools view -s NA21123,NA20752 $i | bcftools norm -m +snps | bcftools view -o debug/vcfcollectionssmall2/$ibase
    bcftools index -t debug/vcfcollectionssmall2/$ibase
done
```

### Recombination maps

Download from [beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip) and change the chromosome names:

```
files=$(ls debug/recombmaps/*.map)
mkdir -p debug/recombmaps2
for i in $files;
do
    ibase=$(basename $i)
    awk '{print("chr"$0)}' $i > debug/recombmaps2/${ibase}
done
```

### DNM

```
cd debug
wget https://static-content.springer.com/esm/art%3A10.1038%2Fnature24018/MediaObjects/41586_2017_BFnature24018_MOESM2_ESM.zip
unzip 41586_2017_BFnature24018_MOESM2_ESM.zip
cd nature24018-s2
tar -xvzf Aging_Oocytes_Supplementary_Table_DNMs.tar.gz
mv decode_DNMs/ ..
cd .. ## back to debug folder
```

### Finally runnning the tool!

```
cargo install --path .
meiosim \
    -r debug/recombmaps2/ \
    -v debug/vcfcollectionssmall2/   \
    -d debug/decode_DNMs/vcfs/Proband-1354.vcf \
    -p NA21123 \
    -P NA20752 \
    --prefix testout \
    --seed 3 -f 5 \
    --genome debug/hg38.genome 
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
