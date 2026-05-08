## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#plots based on the final filtered dataset (f5)

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#load final filtered and annotated snRNA object
sc_IVIg<-qread("IVIG_snRNA_harmony_f5.qs")
Idents(sc_IVIg)<-"celltypes.short"

#=========================================================================================
# nice UMAPS
#=========================================================================================
DimPlot(sc_IVIg,group.by = "celltypes.short.ext",cols = col.ct.ext)
ggsave(filename = "umap annotation major ct ext_leg.jpeg", width=9 , height = 8)

DimPlot(sc_IVIg,group.by = "subcelltypes.short",cols = col.subclusters)
ggsave(filename = "umap annotation subcluster_leg.jpeg", width=15 , height = 8)

#=========================================================================================
# QC plots
#=========================================================================================

dittoBarPlot(sc_IVIg,var = "celltypes.short.ext",group.by = "orig.ident",scale = "count",color.panel = col.ct.ext,retain.factor.levels = T)+coord_flip()
ggsave(filename = "QC plots of final data _ ct count per sample_h.jpeg", width=6 , height = 9)
ggsave(filename = "QC plots of final data _ ct count per sample_h.svg", width=6 , height = 9)

VlnPlot(sc_IVIg,group.by="orig.ident",feature="nFeature_RNA",cols = col.ivig.condi.long,pt.size=0)+NoLegend()
ggsave(filename = "QC plots of final data _  nFeature_RNA.jpeg", width=10 , height = 5)
ggsave(filename = "QC plots of final data _  nFeature_RNA.svg", width=10 , height = 5)

VlnPlot(sc_IVIg,group.by="orig.ident",feature="nCount_RNA",cols = col.ivig.condi.long,pt.size=0)+NoLegend()
ggsave(filename = "QC plots of final data _  nCount_RNA.jpeg", width=10 , height = 5)
ggsave(filename = "QC plots of final data _  nCount_RNA.svg", width=10 , height = 5)

VlnPlot(sc_IVIg,group.by="orig.ident",feature="percent.mt",cols = col.ivig.condi.long,pt.size=0)+NoLegend()
ggsave(filename = "QC plots of final data _  percent.mt.jpeg", width=10 , height = 5)
ggsave(filename = "QC plots of final data _  percent.mt.svg", width=10 , height = 5)

##=========================================================================================
# major cell type markers
#=========================================================================================
#plot of selected classic marker genes per ct in major_extended annoatation
top1 <- c("RYR2","VWF","DCN","ABCC9","MYH11","VCAN","SIGLEC1","FCER1A","KIT","CD3E","KLRF1","NRXN1","MMRN1","NPR3","GPAM")
DotPlot(sc_IVIg, group.by = "celltypes.short.ext",features = top1,dot.scale = 10,assay = "RNA")+ scale_color_gradientn(colors=gradient.col) +coord_flip()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0("top1 markers_dp.jpeg"), width=6.8 , height = 4.5)
ggsave(filename = paste0("top1 markers_dp.svg"), width=6.8 , height = 4.5)

#=========================================================================================
#Plots of gene expression (genes of interest GOI)
#=========================================================================================
#FC receptors
#FcyRIA, FcyRIIA, FcyRIIC, FcyRIIIA (CD16), FcyRIIIB (CD16b), FcyRIIB, FcRn, DC-SIGN
#FCGR1A, FCGR2A, FCGR2C, FCGR3A, FCGR3B, FCGR2B, FCGRT, CD209
FcyR_genes<-c("FCGR1A", "FCGR2A", "FCGR2C", "FCGR3A", 'FCGR3B', "FCGR2B", "FCGRT", "CD209")
sc_IVIg<-ScaleData(sc_IVIg,features = FcyR_genes)

#heatmap + dot
GroupHeatmap(sc_IVIg,features = FcyR_genes,group.by = "celltypes.short",add_dot = T,add_bg = T,
             width = 1.5,height = 1,group_palcolor = list(col.cluster),cell_split_palcolor = list(col.ivig.condi))
ggsave(filename = paste0("Fc receptor expressions in major ct.jpeg") , width=8 , height = 4)
ggsave(filename = paste0("Fc receptor expressions in major ct.svg") , width=8 , height = 4)

