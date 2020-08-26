


cd /home/groups/projects/Sun_RNA_seq/GTEx/Pancreas
ls -lh boot_*_5e+05_long_100/*time* | wc -l
ls -lh boot_*_5e+05_short_100/*time* | wc -l

squeue -u zhabotyn | wc -l
squeue -u zhabotyn | grep ' PD ' | wc -l
squeue -u zhabotyn | grep ' R ' | wc -l
squeue -u zhabotyn | grep ' S ' | wc -l
squeue -u zhabotyn | grep 'p01'



setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/
sbatch -p p01_bat -t 14-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.out --mem=3g --wrap="R360 CMD BATCH '--args 6 246 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer.R /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout"

squeue -u zhabotyn | grep '4811855'


setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/
sbatch -p p01_bat -t 14-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_22_209_371_5e+05_long.out --mem=2g --wrap="R360 CMD BATCH '--args 22 209 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer.R /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_22_209_371_5e+05_long.Rout"

squeue -u zhabotyn | grep '4812048'


setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/
sbatch -p p01_bat -t 14-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.out --mem=2g --wrap="R360 CMD BATCH '--args 8 18 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer.R /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout"

squeue -u zhabotyn | grep '4812149'



setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/
sbatch -p p01_bat -t 14-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.out --mem=2g --wrap="R360 CMD BATCH '--args 4 591 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer.R /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout"

squeue -u zhabotyn | grep '4812263'


tail /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/rout_Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step2_trecase1_22_209_371_5e+05_long.Rout
9664/9664
1922/1922

squeue -u zhabotyn | grep '4812149'
squeue -u zhabotyn | grep '4812263'
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rout_Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rout_Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout

10837/17391
981/2083

tail /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/sstep2_trecase1_22_209_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout

621/9664
1701/1922
6184/17391
339/2083


sbatch -p p01_bat -t 07-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/step5_collectU_5e+05.out --mem=24g --wrap="R CMD BATCH '--args 5e+05' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step5_collectU.R /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/step5_collectU_5e+05.Rout"
Submitted batch job 5045680
squeue -u zhabotyn | grep '5045680'

sbatch -p bat -t 07-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/bout_Brain_Cortex/step5_collectU_5e+05.out --mem=24g --wrap="R CMD BATCH '--args 5e+05' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step5_collectU.R /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/rout_Brain_Cortex/step5_collectU_5e+05.Rout"
Submitted batch job 4995236
squeue -u zhabotyn | grep '4995236'


sbatch -p p01_bat -t 07-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/bout_Whole_Blood/step5_collectU_5e+05.out --mem=24g --wrap="R CMD BATCH '--args 5e+05' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step5_collectU.R /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/rout_Whole_Blood/step5_collectU_5e+05.Rout"
Submitted batch job 4861025
squeue -u zhabotyn | grep '4861025'


sbatch -p p01_bat -t 07-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Pancreas/bout_Pancreas/step5_collectU_5e+05.out --mem=24g --wrap="R CMD BATCH '--args 5e+05' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step5_collectU.R /home/groups/projects/Sun_RNA_seq/GTEx/Pancreas/rout_Pancreas/step5_collectU_5e+05.Rout"
Submitted batch job 4850694
squeue -u zhabotyn | grep '4850694'
tail /home/groups/projects/Sun_RNA_seq/GTEx/Pancreas/rout_Pancreas/step5_collectU_5e+05.Rout
169/194

setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
sbatch -p p01_bat -t 07-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/bout_Whole_Blood/step2_trecase1_6_324_670_5e+05_short.out --mem=4g --wrap="R360 CMD BATCH '--args 6 324 5e+05 short' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1.R /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/rout_Whole_Blood/step2_trecase1_6_324_670_5e+05_short.Rout"
Submitted batch job 4853514
squeue -u zhabotyn | grep '4853514'

tail /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/rout_Whole_Blood/step2_trecase1_6_324_670_5e+05_short.Rout

12408/12409

tail /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/rout_Whole_Blood/step5_collectU_5e+05.Rout
171/171 


tail /home/groups/projects/Sun_RNA_seq/GTEx/Breast_Mammary_Tissue/rout_Breast_Mammary_Tissue/step5_collectU_5e+05.Rout
217/217

tail /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cerebellar_Hemisphere/rout_Brain_Cerebellar_Hemisphere/step5_collectU_5e+05.Rout
90/215

tail /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/rout_Brain_Cortex/step5_collectU_5e+05.Rout
126/211           
tail /home/groups/projects/Sun_RNA_seq/GTEx/Brain_Cortex/step5_collectU_5e+05.Rout
98/221


/home/groups/projects/Sun_RNA_seq/GTEx/Brain_Frontal_Cortex_BA9/rout_Brain_Frontal_Cortex_BA9/step5_collectU_5e+05.Rout
3/210


tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step5_collectU_5e+05.Rout




setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/
sbatch -p bat -t 21-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rerun2_step2_trecase1_8_18_371_5e+05_long.out --mem=4g --wrap="R360 CMD BATCH '--args 8 18 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer2.R /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rerun2_step2_trecase1_8_18_371_5e+05_long.Rout"

setenv LD_LIBRARY_PATH /usr/local/gcc-7.2/lib64
cd /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/
sbatch -p bat -t 21-00:00:00 -o /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rerun2_step2_trecase1_4_591_385_5e+05_long.out --mem=4g --wrap="R360 CMD BATCH '--args 4 591 5e+05 long' /home/groups/projects/Sun_RNA_seq/GTEx/lib/step2_trecase1rer2.R /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rerun2_step2_trecase1_4_591_385_5e+05_long.Rout"

squeue -u zhabotyn | grep '5087702'
squeue -u zhabotyn | grep '5087703'


#rerun1
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout

17391/17391 done
1877/2083   ~88%  12,19

#rerun2
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rerun2_step2_trecase1_8_18_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rerun2_step2_trecase1_4_591_385_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rerun2_step2_trecase1_8_18_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rerun2_step2_trecase1_4_591_385_5e+05_long.Rout

17391/17391 done
2083/2083   done



squeue -u zhabotyn | grep '4812149'
squeue -u zhabotyn | grep '4812263'
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rout_Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/rout_Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/rout_Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout

10837/17391
1194/2083



tail /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/sstep2_trecase1_22_209_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/sstep2_trecase1_22_209_371_5e+05_long.Rout
1922/1922
tail /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout
9664/9664


ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Whole_Blood/step2_trecase1_6_246_670_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/sstep2_trecase1_22_209_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Atrial_Appendage/step2_trecase1_8_18_371_5e+05_long.Rout
ls -lh /home/groups/projects/Sun_RNA_seq/GTEx/Heart_Left_Ventricle/step2_trecase1_4_591_385_5e+05_long.Rout
