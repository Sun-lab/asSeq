cut -f 3 Mus_musculus.NCBIM37.67.gtf | sort | uniq -c > "step0_check_data.log"

cut -f 1 Mus_musculus.NCBIM37.67.gtf | sort | uniq -c >> "step0_check_data.log"

# cut -f 1 Mus_musculus.NCBIM37.67.updated.gtf | sort | uniq -c
