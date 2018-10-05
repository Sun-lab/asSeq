# ----------------------------------------------------------------------
# check if the snps labeled as "strand issue"
# can be fixed by flipping them with plink
# ----------------------------------------------------------------------
# the location of shapeit and reference panel
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
plinkloc = "plink"

# the location of chromoson level data and precheck in step1
proj = "../../data/chr"
gwas_alignments = "../../data/chr/precheck_log"

checki <- matrix(NA, nrow = 22, ncol = 2)
rownames(checki) = sprintf("chr%s", 1:22)

#as with each shapeit we do it on chromosome level
all1000G = "1000GP_Phase3" #common reference prefix
for (i in seq(22)) {
  #i = 22;
  #flip the snps of interest
  con = file(sprintf("%s/gwas.alignments_chr%d.snp.strand", gwas_alignments, i))
  lns = readLines(con)
  lns = lns[substr(lns, 1, 6) == "Strand"]
  strand = matrix(unlist(strsplit(lns, split = "\t")), nrow = 11)[4, ]
  close(con)
  write.table(
    unique(strand),
    "tmp_strand.txt",
    row.names = F,
    col.names = F,
    quote = F
  )
  file.copy(from = sprintf("%s/geno.chr%s.map", proj, i),
            to = "tmp.map")
  file.copy(from = sprintf("%s/geno.chr%s.ped", proj, i),
            to = "tmp.ped")
  system(sprintf("%s --file tmp --flip tmp_strand.txt  --recode", plinkloc))
  cat("plink done ")
  
  con = file(sprintf("tmp.ped", proj, i))
  tmp = readLines(con)
  close(con)
  cat("\t read preconv ")
  
  con = file(sprintf("plink.ped", proj, i))
  tmp0 = readLines(con)
  close(con)
  cat("\t read postconv: ")
  #check that the plink really flipped those snps
  #we see that plink does flip those strands
  message(!all(unlist(tmp) == unlist(tmp0)))
  #recheck using shapeit on flipped data
  pref = sprintf("%s -check", shapeloc)
  pedi = sprintf("--input-ped plink")
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",
                 refinfo,
                 i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz", refinfo, all1000G, i)
  legr = sprintf("%s/%s_chr%s.legend.gz", refinfo, all1000G, i)
  samr = sprintf("%s/%s.sample", refinfo, all1000G)
  out = sprintf("--output-log tmp")
  com = sprintf("%s %s %s %s %s %s %s",
                pref, pedi, mapr, inpr, legr, samr, out)
  message(com)
  system(com)
  
  
  con = file(sprintf("tmp.snp.strand", proj, i))
  lns2 = readLines(con)
  lns2 = lns2[substr(lns2, 1, 6) == "Strand"]
  strand2 = matrix(unlist(strsplit(lns2, split = "\t")), nrow = 11)[4, ]
  close(con)
  
  #count how many snps are in the one exclusion set, but not in the other
  checki[i, 1] = length(setdiff(unique(strand), unique(strand2)))
  checki[i, 2] = length(setdiff(unique(strand2), unique(strand)))
  message("processed ", i, "'th file")
}

#remove temporary files
system("rm -f tmp*")
system("rm -f plink*")

checki

sessionInfo()
q("no")
