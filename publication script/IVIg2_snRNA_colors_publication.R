## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#hand picked colors for the project
library(RColorBrewer)
library(rcartocolor)

#------function to get the standart colours of seurat-----------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#12 colorblind safe colours from rcartocolor
#display_carto_all(colorblind_friendly = TRUE)
safe_c=carto_pal(12,"Safe")
safe_c2=safe_c[c(1:8,12,10,11,9)] #change order so the two darker ones are not next to each other


#-----------colours for the plots
# project condition colours
col.ivig.condi<-c("#B2DF8A","#33A02C","#CAB2D6","#6A3D9A")
col.ivig.condi.long<-c(rep(col.ivig.condi[1],10),rep(col.ivig.condi[2],11),rep(col.ivig.condi[3],9),rep(col.ivig.condi[4],9))

#gradient from Placebo to IVIg
gradient.col.condi = c("#33A02C","white","#6A3D9A")

#simple 24 colors vector, one color per sample
sample.cols<-gg_color_hue(42)
#cluster colors
col.cluster<-safe_c2[c(1,3,2,4:12)]
names(col.cluster)<-c("CM","EC","FIB","Mural","Myeloid","Lymphoid", "EndoEC","LYMEC","NEURO","ADI" )

#subcluster colors
col.CM<-colorRampPalette(c("black",safe_c2[1], "white"))
col.CM<-col.CM(8)[2:7]
col.CM<-col.CM[-6]
col.EC<-colorRampPalette(c("black",safe_c2[3], "white"))
col.EC<-col.EC(6)[2:5]
col.EC<-c(col.EC,"#44AA99","#999933")
col.FIB<-colorRampPalette(c("black",safe_c2[2], "white"))
col.FIB<-col.FIB(8)[2:7]
col.FIB<-col.FIB[c(6,1,5,3,4,2)]

#myeloid cols
col.Myel<- colorRampPalette(c("white","#0000FF","grey", safe_c2[5], "#88CCEE","black"))(16)
col.Myel<-col.Myel[c(2:6,8:15)]
col.Myel<-col.Myel[c(3,2,1,4:13)]
cell.levels.old = c("MP_LYVE1hi_MHCIIlow","MP_LYVE1hi_MHCIIint","MP_LYVE1low_MHCIIhi","cDC","MO_CD16","MO_VCAN","MP_FOLR2","MP_TREM2","MP_ISG","MP_NFKB","MAST","MP_proliferating1","MP_proliferating2")
names(col.Myel)<-cell.levels.old
cell.levels = c("cDC","MO_CD16","MO_VCAN","MP_LYVE1hi_MHCIIlow","MP_LYVE1hi_MHCIIint","MP_LYVE1low_MHCIIhi","MP_FOLR2","MP_TREM2","MP_ISG","MP_NFKB","MP_proliferating1","MP_proliferating2","MAST")
col.Myel<-col.Myel[cell.levels]
names(col.Myel)<-NULL

col.Lymph<-colorRampPalette(c("black",safe_c2[6], "white"))
col.Lymph<-col.Lymph(12)[2:11]
col.Lymph<-col.Lymph[c(1,6,3,4,5,9,7,8,2,10)]
col.Mural<-colorRampPalette(c("black",safe_c2[4], "white"))
col.Mural<-col.Mural(11)[2:10]
col.Mural<-col.Mural[c(5,1,9,4,2,8,6,3)]
col.NEURO<-colorRampPalette(c("black",safe_c[12], "white"))
col.NEURO<-col.NEURO(5)[2:4]
col.ADI<-safe_c[10]

col.subclusters<-c(col.CM,col.EC[1:4],col.FIB,col.Mural,col.Myel,col.Lymph,col.EC[5:6],col.NEURO,col.ADI)

col.patient<-list(rep(col.ivig.condi[1],8),rep(col.ivig.condi[2],11),
               rep(col.ivig.condi[3],9),rep(col.ivig.condi[4],9))

col.ct.ext<-c(col.cluster[1:3],
              colorRampPalette(c("black","#117733", "white"))(4)[c(2,3)],
              colorRampPalette(c("black","#0000FF", "white"))(6)[c(2:5)],
              colorRampPalette(c("black",safe_c2[6], "white"))(4)[c(2,3)],
              col.cluster[c(9,8,7,10)])
names(col.ct.ext)<-NULL

#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))



