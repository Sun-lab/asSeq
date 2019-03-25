library("GenomicFeatures")

proj_input = "/fh/fast/sun_w/licai/_tumor_eQTL"
gtfIn   = sprintf("%s/gencode.v28.annotation.gtf", proj_input)  # from GENCODE

txdb2 = makeTxDbFromGFF(file=gtfIn, format="gtf",
                        # dataSource=paste(path, file, sep=""),
                        organism="Homo sapiens")

saveDb(txdb2, file=sprintf("%s/Homo_sapiens_gencode_v28_GRCh38.sqlite", proj_input))

sessionInfo()
q(save="no")
