#!/bin/bash
# Example invocation of the Regulatory Locus Intersection (RELI) tool using
# sample data in the 'example' subdirectory (European ancestry)

RELIBIN=${RELI_BINARY:-../RELI}
DATADIR=${RELI_DATA_DIR:-../data}
OUTPUTDIR=${RELI_OUTPUT_DIR:-../output}

OPTIONS=(
    -index "$DATADIR/ChIPseq.index"
    -data "$DATADIR/ChIP-seq"
    -target hg19_0302
    -build "$DATADIR/GenomeBuild/hg19.txt"
    -out "$OUTPUTDIR"
    -species hg19
    -null OpenChrom
    -rep 2000
)

# invoke RELI; double quotes are essential here in case any of the options
# above contain embedded spaces
"$RELIBIN" "${OPTIONS[@]}"

