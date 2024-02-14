

build:
    cargo build --release && cargo build

clean:
    rm -rf testout*

run prefix:
    #!/usr/bin/env sh
    echo 'Running {{prefix}}…'
    time ./target/release/meiosim \
        -r debug/recombmaps2/ \
        -v debug/vcfcollectionssmall2/      \
        -d debug/decode_DNMs/vcfs/Proband-1354.vcf     \
        -p NA21123     \
        -P NA20752  \
        --seed 10 \
        -f 5     \
        --genome debug/hg38.genome     \
        --dwgsim \
        --prefix {{prefix}}

lines:
    wc -l testout*/*
