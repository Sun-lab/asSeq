The basic pipeline includes data preparation steps:

1 step0_submit_trim.R - trimming suspicious allele-specific counts and total counts

2 step1_submit_preprocSNP.R - creating separate files for each gene

main eQTL analysis and p-value estimate:

3 step2_submit_trecaseA.R

4 step4_submit_MatrixEQTL.R

and further analysis steps including principal component analysis and dynamic eqtl findigs.

Once the specification file is set up, each script should be able to be run authomatically (with corresponding note given to the script to distinguish between long and short model, etc)

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt long' step0_submit_trim.R MS_step0_submit_trim_long.Rout

#to save time short model fit relies on long model being run first, so this step should be done after the previous is completed

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt short' step0_submit_trim.R MS_step0_submit_trim_short.Rout

#wait until the previous step completes

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt' step1_submit_preprocSNP.R MS_step1_submit_preprocSNP.Rout

#wait until the previous step completes

#fitting steps can be run at once, if you have generous enough scheduler 

#(often schedulers have limit on the number of maximum number of jobs that can be submitted)

#fitting TReCASE model for each gene, long model:

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt long 5e5' step2_submit_trecaseA.R MS_step2_submit_trecaseA_long.Rout

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt short 5e5' step2_submit_trecaseA.R MS_step2_submit_trecaseA_short.Rout

#estimating permuted p-values

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt long 5e5' step4_submit_MatrixEQTL.R MS_step4_submit_MatrixEQTL_long.Rout

R CMD BATCH  '--args specifications_Muscle_Skeletal.txt short 5e5' step4_submit_MatrixEQTL.R MS_step4_submit_MatrixEQTL_short.Rout

