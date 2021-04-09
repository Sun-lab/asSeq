
library(ggplot2)
library(ggrepel)

theme_set(theme_classic())

# -------------------------------------------------------------------------------------
# set up colors
# -------------------------------------------------------------------------------------

tiscol = c(
rgb(170, 170, 255, maxColorValue = 255),#MS , #AAAAFF - light blue
rgb(255,0,187, maxColorValue = 255), #Whole Blood - rgb() - #FF00BB - dark pink
rgb(119, 119, 255, maxColorValue = 255),#SSELl - rgb() - #7777FF - light blue
rgb(255, 0, 0, maxColorValue = 255),#Ar. Tib -   - #FF0000
rgb(255, 102, 0, maxColorValue = 255),#Ad. Sub  -     - #FF6600
rgb(0,102,0, maxColorValue = 255),#Thyroid - rgb() - #006600 - dark green
rgb(255, 215, 0, maxColorValue = 255), #NT- - #FFD700 - light yellow
rgb(0, 0, 255, maxColorValue = 255),#SNSES - rgb() - #0000FF - dark blue
rgb(153, 255, 0, maxColorValue = 255), #Lung - #99FF00 - light green
rgb(85, 34, 0, maxColorValue = 255),#EMuc -  - #552200 - dark orange
rgb(170, 238, 255, maxColorValue = 255),#CCf      -  - #AAEEFF
rgb(255, 170, 0, maxColorValue = 255),#Ad. Visc -     - #FFAA00
rgb(136, 153, 136, maxColorValue = 255),#EMus -  - #889988 - light green
rgb(51,204,204, maxColorValue = 255),#Brst.MT  -  - #33CCCC
rgb(255, 85, 85, maxColorValue = 255),#Ar. Aort -     - #FF5555
rgb(102, 0, 153, maxColorValue = 255),#HLV - #660099 - dark violet
rgb(153, 0, 255, maxColorValue = 255),#HAA -  - #9900FF - dark violet
rgb(204, 153, 85, maxColorValue = 255),#CT       -  - #CC9955
rgb(136, 115, 85, maxColorValue = 255),#EGJ -  - #887355 - dark orange
rgb(255, 221, 153, maxColorValue = 255),#Stomach - rgb() - #FFDD99 - light orange
rgb(170,170,170, maxColorValue = 255), #Testis - rgb() - #AAAAAA - light red
rgb(153, 85, 34, maxColorValue = 255),#Panc - rgb() - #995522 - dark orange
rgb(255, 170, 153, maxColorValue = 255),#Ar. Coro - 	- #FFAA99
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255)#Br       -     - #EEEE00
)
tiss = c(
"Muscle_Skeletal",
"Whole_Blood",
"Skin_Sun_Exposed_Lower_leg",
"Artery_Tibial",
"Adipose_Subcutaneous",
"Thyroid",
"Nerve_Tibial",
"Skin_Not_Sun_Exposed_Suprapubic",
"Lung",
"Esophagus_Mucosa",
"Cells_Cultured_fibroblasts",
"Adipose_Visceral_Omentum",
"Esophagus_Muscularis",
"Breast_Mammary_Tissue",
"Artery_Aorta",
"Heart_Left_Ventricle",
"Heart_Atrial_Appendage",
"Colon_Transverse",
"Esophagus_Gastroesophageal_Junction",
"Stomach",
"Testis",
"Pancreas",
"Artery_Coronary",
"Brain_Cortex",
"Brain_Nucleus_accumbens_basal_ganglia",
"Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere","Brain_Frontal_Cortex_BA9"
)

abbs = c("MS", "WB", "SSELl", "AT", "AS", "Th", "NT", "SNSES", "L", "EMuc", "CCf", 
"AVO", "EMus", "BMT", "AA", "HLV", "HAA", "CT", "EGJ", "S", "Te", "P", "AC", "BC",
"BNabg", "BCbg", "BCH", "BFC")

nsam = c(704, 670, 603, 581, 579, 573, 532, 514, 513, 495, 482, 467, 463, 396,
386, 385, 371, 366, 329, 323, 321, 303, 212, 205, 202, 194, 175, 175)

tiscol = tiscol[order(tiss)]

# -------------------------------------------------------------------------------------
# plot for one variable of dynamic eQTL
# -------------------------------------------------------------------------------------

plot_n_dyn_eQTL <- function(dv){
  # dv = "age"
  
  dat = read.csv(sprintf("dyn_cand/summary_significant_%s.csv", dv))
  dim(dat)
  dat[1:2,]
  
  setequal(dat$X, tiss)
  
  ti = dv
  if(dv != "age"){ti = toupper(dv) }
  
  g1 = ggplot(dat, aes(x=S0.10+1, y=L0.10, color=X)) +
    geom_point() + scale_color_manual(values=tiscol) + 
    theme(legend.title=element_blank()) + 
    labs(x = "# of dyn_eQTLs from short model", title = ti, 
         y = "# of dyn_eQTLs from long model") + 
    scale_x_continuous(trans='log10', breaks=c(1, 11, 21, 101), 
                       labels=c("0", "10", "20", "100"))
  
  pdf(sprintf("n_dyn_eGenes_%s_with_legend.pdf", dv), width=10, height=3.3)
  print(g1)
  dev.off()
  
  dat$tissue = rep("", nrow(dat))
  if(dv == "tp53"){
    w2label    = which(dat$S0.10 > 200)
  }else{
    w2label    = which(dat$S0.10 > 100)
  }
  dat$tissue[w2label] = dat$X[w2label]
  
  g2 = g1 + theme(legend.position = "none") + 
    geom_text_repel(aes(label = tissue), data = dat, size = 3.5, force=2)
  
  pdf(sprintf("n_dyn_eGenes_%s.pdf", dv), width=2.8, height=2.8)
  print(g2)
  dev.off()
}

plot_n_dyn_eQTL("age")
plot_n_dyn_eQTL("ctcf")
plot_n_dyn_eQTL("tp53")


sessionInfo()
q(save="no")
