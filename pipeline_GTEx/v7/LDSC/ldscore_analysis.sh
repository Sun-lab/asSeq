#!/bin/bash
ct=$1
echo 'cell type:' ${ct}
cd ../ldsc/
module load anaconda2
source activate ldsc
ml bedtools
python ldsc.py \
	--h2 ../data/clozuk_pgc2.meta.sumstats_with_rsid.gz \
	--w-ld-chr weights_hm3_no_hla/weights. \
	--ref-ld-chr ../data_eQTL/${ct}.,../data_eQTL/control_${ct}. \
	--frqfile-chr 1000G.mac5eur. \
	--out ../output_eQTL/${ct} \
	--print-coefficients
source deactivate
