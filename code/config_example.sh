#!/usr/bin/env bash

RESULTS=results
PENNCNV_RESULTS=$RESULTS/penncnv
DNACOPY_RESULTS=$RESULTS/dnacopy
VICE_RESULTS=$RESULTS/vanilla_ice

if [[ ! -d $PENNCNV_RESULTS ]]; then
        mkdir -p $PENNCNV_RESULTS
fi
if [[ ! -d $DNACOPY_RESULTS ]]; then
        mkdir -p $DNACOPY_RESULTS
fi
if [[ ! -d $VICE_RESULTS ]]; then
        mkdir -p $VICE_RESULTS
fi
if [[ ! -d logs ]]; then
        mkdir logs
fi


# For testing we will use the example data that comes with PennCNV

penncnv_folder=`which detect_cnv.pl`
$penncnv_folder=`sed 's/\/detect_cnv.pl//'`

EXAMPLE_DATA=(
$penncnv_folder/example/mother.txt
$penncnv_folder/example/father.txt
$penncnv_folder/example/offspring.txt
$penncnv_folder/example/example.pfb
$penncnv_folder/example/example.hmm
$penncnv_folder/example/GenomeWideSNP_6.hg18.cnv_defs
)

for x in ${EXAMPLE_DATA[@]}
do
	if [ ! -f $x ]; then
    	echo "Cannot locate the example data from the PennCNV folder"
    	echo "Please ensure PennCNV is installed and has been added to your \$PATH"
    	exit 1
	fi
done

SAMPLES=(
$penncnv_folder/example/mother.txt
$penncnv_folder/example/father.txt
$penncnv_folder/example/offspring.txt
)



