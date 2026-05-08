## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#initial assessment of data QC and subsequent filtering for high quality nuclei transcriptoms

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")
cellranger_out_path <- "./cellranger_output"

#=========================================================================================
#load sample meta data
#=========================================================================================
meta.data <- as.data.frame(read_xlsx("sample_meta_data_publication_reordered.xlsx"))

#=========================================================================================
#loop though all cell ranger out folders and perform basic processing per sample
#=========================================================================================
#load folder list
sample_folder_list<-list.files(cellranger_out_path)

#check list
sample_folder_list %in% meta.data$study_sample_ID
meta.data$study_sample_ID %in% sample_folder_list

#adjust sample folder list order (important for reproducibility)
sample_folder_list<-meta.data$study_sample_ID

#----------loop through all individual samples for QC-----------
#QC check for each sample
#1. run clustering of the individual sample as this is required for SoupX and scDblFinder 
#2. SoupX to find ambient RNA contamination
#3. scDblFinder to unbiasedly check for doublets
#4. plot results for each sample separate (combined PDF)
  #bad sample: "IVIG_3_POST"
  which("IVIG_3_POST"==sample_folder_list)
  
  filtered.samples.list=NULL
  plots.of.filtering=NULL
  for (sample_ID in sample_folder_list[-12]) {
    print(sample_ID)
    #load data and create seurat object
    path<-paste0(cellranger_out_path,sample_ID,"/outs/filtered_feature_bc_matrix/")
    if (dir.exists(path)) {
      sc.data <- Read10X(data.dir = path)
      sc <- CreateSeuratObject(counts=sc.data, project=sample_ID,min.cells = 3, min.features = 200)
      
      #assess quality by function
      sample_path <-paste0(cellranger_out_path,sample_ID,"/outs")
      qc.result <- assess_quality(sc.object = sc, res =0.1,
                                   runSoupX=T,
                                   sample_path=sample_path,
                                   upperFeatureCutoff = 6000,lowerFeatureCutoff=500,mtCutoff =15)
      filtered.samples.list[[sample_ID]]=qc.result[[1]]
      plots.of.filtering[[sample_ID]]=qc.result[[2]]

      #plot one pdf per sample
      pdf(file = paste0(sample_ID,"_quality check.pdf"),paper = "a4",height = 23,width = 33)
      title <- ggdraw() + draw_label(sample_ID,fontface = 'bold',x = 0,hjust = 0)+theme(plot.margin = margin(0, 0, 0, 7))
      print(plot_grid(title, qc.result[[2]],ncol = 1,rel_heights = c(0.03, 1)))
      dev.off()

    }#end of if path
  }

  #save all assessed seurat objects
  qsave(filtered.samples.list,file = "list of individual samples.qs")
  #filtered.samples.list_p<-qread("./list of individual samples.qs")
  
  #plot QC dimplot overview
  qc.overview <-NULL
  for (i in 1:length(filtered.samples.list)) {
    qc.overview[[i]]<- DimPlot(filtered.samples.list[[i]],group.by = "QC",pt.size = 0.5,cols = rev(dittoColors(1)[1:5]),
                               order = rev(c("Pass","Low_nFeature","High_nFeature","High_MT","Doublet")))+
      theme(title = element_blank(),axis.text = element_blank(),text=element_text(size=8),
            legend.key.height= unit(0.1, 'cm'),legend.key.width= unit(0.1, 'cm'),
            legend.position = c(0.6,0.9))
  }
  ggarrange(plotlist = qc.overview,common.legend = T,legend = "bottom")
  ggsave(filename = "sample QC overview.jpeg", width=35 , height = 30,limitsize = F)
  ggsave(filename = "sample QC overview.svg", width=35 , height = 30,limitsize = F)
  
#=========================================================================================
# merge sample, check QC and run harmony. filter low quality nuc and doublets
#=========================================================================================
  
#-------merge the filtered samples back into one object
  Samples.combined<-merge(filtered.samples.list[[1]],filtered.samples.list[-1])
  rm(filtered.samples.list)
  gc()
  
  #convert to an v3 seurat object (v5 causes many problems with other functions)
  Samples.combined[["RNA"]] <- as(object = Samples.combined[["RNA"]], Class = "Assay")
  
#--------run harmony batch correction
  Samples.combined <- NormalizeData(Samples.combined, verbose = T)
  Samples.combined <- FindVariableFeatures(Samples.combined, verbose = T)
  Samples.combined <- ScaleData(Samples.combined, verbose = T)
  Samples.combined <- RunPCA(Samples.combined, verbose = T)
  ElbowPlot(Samples.combined,ndims = 50)
  Samples.combined <- RunUMAP(Samples.combined,dims = 1:20, verbose = T)

  #plot batch effect
  DimPlot(Samples.combined,group.by = "orig.ident",cols = sample.cols)
  ggsave(filename = "batch effect raw.jpeg", width=15 , height = 10)
  DimPlot(Samples.combined,group.by = "orig.ident",cols = sample.cols)+NoLegend()
  ggsave(filename = "batch effect raw_nolegend.jpeg", width=10 , height = 10)
  
  #run harmony
  Samples.combined <- Samples.combined %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  ElbowPlot(Samples.combined,reduction = "harmony",ndims = 50)
  #process based on harmony result
  Samples.combined <- Samples.combined %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    identity()
  
  #high res
  Samples.combined <- Samples.combined %>%  FindClusters(resolution = 1)
  DimPlot(Samples.combined,label = T,label.size = 10)+NoLegend()
  ggsave(filename = "harmony umap res 1.jpeg", width=10 , height = 10)
  
  dittoBarPlot(Samples.combined,var = "QC",group.by = "seurat_clusters",retain.factor.levels = T)
  ggsave(filename = "per seuratcluster QC.jpeg", width=10 , height = 5)
  ggsave(filename = "per seuratcluster QC.svg", width=10 , height = 5)
  
  #common major cell type markers
  Samples.combined<-ScaleData(Samples.combined,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  FeaturePlotCombined(Samples.combined,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  ggsave(filename = "after QC and harmony.jpeg", width=20 , height = 20)
  
  #plot MT-count unfiltered
  FeaturePlotCombined(Samples.combined,features = c("nFeature_RNA"))
  ggsave(filename = "harmony raw nFeature_RNA.jpeg", width=10 , height = 10)
  FeaturePlotCombined(Samples.combined,features = c("percent.mt"))
  ggsave(filename = "harmony raw percent.mt.jpeg", width=10 , height = 10)

  VlnPlot(Samples.combined,group.by="seurat_clusters",feature="nFeature_RNA",pt.size=0)+NoLegend()
  ggsave(filename = "nFeature_RNA res 1.jpeg", width=10 , height = 5)
  
  VlnPlot(Samples.combined,group.by="seurat_clusters",feature="percent.mt",cols = sample.cols,pt.size=0)+NoLegend()
  ggsave(filename = "percent.mt res1.jpeg", width=10 , height = 5)
  
  #MT-counts is unexpected very low, mostly <5%
  #therefore cutoff set from 15% to 5%
  Samples.combined[['QC']] <- ifelse(Samples.combined$percent.mt > 5,'High_MT',Samples.combined$QC)

  #save merged raw
  qsave(Samples.combined,file = "merged raw _harmony.qs")

#------------find and remove low quality clusters
#low nfeature
  #9, 24, 34

#high MT
  #8, 9, 11, 20, 26, 31 (14, but also has high DCN)

#high doublet score
  #13, 15, 22, 25, 26, 28, 29, 33

#8 looks like CM with EC contamination and high MT    

#low quality cluster
  Idents(Samples.combined)="seurat_clusters"
   cells_rm<-WhichCells(Samples.combined,idents = c(9, 24, 34, #low nfeature
                                                   8, 11, 20, 26, 31, #high MT
                                                   13, 15, 22, 25, 26, 28, 29, 33))#high doublet score
  #also remove all nuc's with more then 5% MT count
  cells_rm<-unique(c(cells_rm,WhichCells(Samples.combined, expression = percent.mt > 5)))
  
  #save with cells removed to be able to highligh them later
  saveRDS(cells_rm,"cells of clusters filtered due to low nfeature higMT and high doublet score.RDS")
  
  #rm
  Samples.combined<-subset(Samples.combined,cells= cells_rm,invert=T)  
  
#----------run harmony pipeline again on filtered
  Samples.combined <- FindVariableFeatures(Samples.combined, verbose = T)
  Samples.combined <- ScaleData(Samples.combined, verbose = T)
  Samples.combined <- RunPCA(Samples.combined, verbose = T)
  Samples.combined <- Samples.combined %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  ElbowPlot(Samples.combined,reduction = "harmony",ndims = 50)
  #process based on harmony result
  Samples.combined <- Samples.combined %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    identity()
  Samples.combined <- Samples.combined %>%  FindClusters(resolution = 0.1)
  
  DimPlot(Samples.combined,label = T,label.size = 10)+NoLegend()
  ggsave(filename = "after f1 res0.1  .jpeg", width=10 , height = 10)

  #common major celltype markers
  Samples.combined<-ScaleData(Samples.combined,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  FeaturePlotCombined(Samples.combined,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  ggsave(filename = "after f1 marker genes.jpeg", width=20 , height = 15)
  
#=========================================================================================
# preliminary annotation of major celltypes for subclustering an further filtering
#=========================================================================================

#------------------------------annotation of major celltype clusters
  Idents(Samples.combined)="seurat_clusters"
  types <- list("0"="CM","1"="FIB","2"="EC","3"="Mural","4"="EndoEC","5"="Myeloid","6"="Lymphoid","7"="NEURO","8"="LYMEC","9"="ADI")
  #Rename identities
  Samples.combined <- RenameIdents(Samples.combined, types)
  Samples.combined$celltypes.short <- Idents(Samples.combined)
  cell.levels = c("CM","EC","FIB","Mural","Myeloid","Lymphoid","EndoEC","LYMEC","NEURO","ADI")
  Samples.combined$celltypes.short = factor(Samples.combined$celltypes.short,levels = cell.levels)
  Idents(Samples.combined)="celltypes.short"

  #plot celltypes in the umap
  DimPlot(Samples.combined,label = T,label.size = 10)+NoLegend()
  ggsave(filename = "after f1 major cell types annotation.jpeg", width=10 , height = 10)

  #save after first filter
  qsave(Samples.combined,file = "IVIG_snRNA_harmony_f1.qs")
  # #load
  # Samples.combined<-qread(file = "IVIG_snRNA_harmony_f1.qs")
  
#----------------loop through major cell types for subclustering and subsequent identification of doublet and low quality subclusters
#often low quality nuclei and doublets can only be identified on subcluster level, therefore each major celltype was subclustered and checked
    
#loop through all cell types
for (celltype in levels(Samples.combined$celltypes.short)) {
  #create out dir
  dir.create(paste0("./",celltype))
  outfld = paste0("./",celltype,"/")
  
  #get subtype and subcluster incl. reintegration by harmony
  Samples.combined.sub<-subset(Samples.combined,idents = celltype)
  Samples.combined.sub <- FindVariableFeatures(Samples.combined.sub, verbose = FALSE)
  Samples.combined.sub <- ScaleData(Samples.combined.sub, verbose = T)
  Samples.combined.sub <- RunPCA(Samples.combined.sub, verbose = T)
  Samples.combined.sub <- Samples.combined.sub %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  Samples.combined.sub <- Samples.combined.sub %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    identity()
  Samples.combined.sub <- Samples.combined.sub %>%  FindClusters(resolution = 0.3)
  
  #save subclustering for separate filtering 
  qsave(Samples.combined.sub,file = paste0(outfld,"subclustering of ",celltype,".qs"))
  
  #plot clusters
  DimPlot(Samples.combined.sub,label = T,label.size = 10,pt.size = 2)+NoLegend()
  ggsave(filename = paste0(outfld,"subclustering of ",celltype,"1.jpeg"), width=5 , height = 5)
  DimPlot(Samples.combined.sub,label = T,label.size = 10)+NoLegend()
  ggsave(filename = paste0(outfld,"subclustering of ",celltype,"2.jpeg"), width=10 , height = 10)
  
  #plot classic markers
  FeaturePlotCombined(Samples.combined.sub,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  ggsave(filename = paste0(outfld,"marker genes ",celltype,".jpeg") , width=20 , height = 20)
  
  #get cluster marker
  all.markers = FindAllMarkers(Samples.combined.sub, min.pct = 0.3)
  saveRDS(all.markers,file = paste0(outfld,celltype,"all.markers.RDS"))
  top10 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene))
  DotPlot(Samples.combined.sub, features = rev(top10),dot.scale = 10,assay = "RNA")+ coord_flip()
  ggsave(filename = paste0(outfld,celltype,"top10 markers.jpeg"), width=10 , height = 30)
  
  FeaturePlotCombined(Samples.combined.sub,features = c("percent.mt","nFeature_RNA"))
  ggsave(filename = paste0(outfld,"mt percent and nFeature ",celltype,".jpeg") , width=10 , height = 5)
  
}

#----------------remove cluster with doublets from the subclustering
#collect all cells that will be removed to rm them from the full merged dataset
cells.to.rm = NULL
#CM
  Samples.combined.CM <- qread("./QC/subclustering/CM/subclustering of CM.qs")
  cells.to.rm<-WhichCells(Samples.combined.CM,idents= c(3))
#Mural
  Samples.combined.Mural <- qread("./QC/subclustering/Mural/subclustering of Mural.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.Mural,idents= c(5,6)))
#EC
  Samples.combined.EC <- qread("./QC/subclustering/EC/subclustering of EC.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.EC,idents= c(6)))
#FIB
  Samples.combined.FIB <- qread("./QC/subclustering/FIB/subclustering of FIB.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.FIB,idents= c(5,6)))
#Myeloid 
  Samples.combined.Myeloid <- qread("./QC/subclustering/Myeloid/subclustering of Myeloid.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.Myeloid,idents= c(3,5,10,11)))
#LYMEC
  Samples.combined.LYMEC <- qread("./QC/subclustering/LYMEC/subclustering of LYMEC.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.LYMEC,idents= c(1,6,4,3,7)))
#NEURO
  Samples.combined.NEURO <- qread("./QC/subclustering/NEURO/subclustering of NEURO.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.NEURO,idents= c(3,4,5,1,9,6,8)))
#EndoEC
  Samples.combined.EndoEC <- qread("./QC/subclustering/EndoEC/subclustering of EndoEC.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.EndoEC,idents= c(4,6)))
#Lymphoid
  Samples.combined.Lymphoid <- qread("./QC/subclustering/Lymphoid/subclustering of Lymphoid.qs")
  cells.to.rm<-c(cells.to.rm,WhichCells(Samples.combined.Lymphoid,idents= c(5,4)))
  
#save list of cells to remove
  saveRDS(cells.to.rm,"subcluster based filtering cells to remove.RDS")
  
#---------------remove doublets and low quality nuc identified in subclustering
  #load
  cells.to.rm<-readRDS("subcluster based filtering cells to remove.RDS")
  #rm
  Samples.combined<-subset(Samples.combined,cells= cells.to.rm,invert=T)  
  
  #plot celltypes in the umap
  DimPlot(Samples.combined,label = T,label.size = 10)+NoLegend()
  ggsave(filename = "after f2 major cell types annotation.jpeg", width=10 , height = 10)
  
#=========================================================================================
# add sample level meta data
#=========================================================================================
#sample meta data
  meta.data
  meta.data$condition<-paste0(meta.data$Study_group,"_",meta.data$time_point)
  meta.data$patient_ID<-meta.data$study_sample_ID
  meta.data$study_sample_ID<-NULL
  #get metadata of snRNA
  Samples.combined$patient_ID<-Samples.combined$orig.ident
  snRNA_meta.data<-Samples.combined@meta.data
  snRNA_meta.data<-left_join(snRNA_meta.data,meta.data,"patient_ID")
  snRNA_meta.data$patient_ID<-as.character(snRNA_meta.data$patient_ID)
  snRNA_meta.data$patient_ID<-gsub("_PRE","",snRNA_meta.data$patient_ID)
  snRNA_meta.data$patient_ID<-gsub("_POST","",snRNA_meta.data$patient_ID)
  rownames(snRNA_meta.data)<-rownames(Samples.combined@meta.data)
  #add metadata
  Samples.combined@meta.data<-snRNA_meta.data
  
  #save after adding meta data
  qsave(Samples.combined,file = "IVIG_snRNA_harmony_f2.qs")

#=========================================================================================
# nFeautre coverage and sex specific gene expression clearly identified problematic samples
#=========================================================================================
  #load
  Samples.combined<-qread("IVIG_snRNA_harmony_f2.qs")
  
  #remove IVIG_9_PRE as it is very clearly a bad sample (very low nFeature count)
  Idents(Samples.combined)<-"orig.ident"
  Samples.combined<-subset(Samples.combined,ident="IVIG_9_PRE",invert=T)
  
  #remove Placebo_9_PRE
  #this sample should be male according to meta data, but has female gene expression signature
  #Placebo_9_PRE (corresponding pair) is male as it should be
  #a sample mix-up is likely, therefore this sample needs to be removed
  Samples.combined<-subset(Samples.combined,ident="Placebo_9_PRE",invert=T)    
  
  #add levels to conditions and sample
  Samples.combined$condition<-factor(Samples.combined$condition,levels = c("PLACEBO_PRE","PLACEBO_POST","IVIG_PRE" ,"IVIG_POST"))
  Samples.combined$time_point<-factor(Samples.combined$time_point,levels = c("PRE","POST"))
  Samples.combined$Study_group<-factor(Samples.combined$Study_group,levels = c("PLACEBO","IVIG"))
  
  meta.data<-Samples.combined@meta.data
  meta.data<-arrange(meta.data,condition)
  unique(meta.data$orig.ident)
  Samples.combined$orig.ident<-factor(Samples.combined$orig.ident,levels = unique(meta.data$orig.ident))
  Samples.combined$patient_ID<-factor(Samples.combined$patient_ID,levels = unique(meta.data$patient_ID))
  
  #save after filter
  qsave(Samples.combined,file = "IVIG_snRNA_harmony_f3.qs")
  
  #plot celltypes in the umap
  DimPlot(Samples.combined,label = T,label.size = 10,group.by = "celltypes.short",cols = safe_c2)+NoLegend()
  ggsave(filename = "after f3 major cell types annotation.jpeg", width=10 , height = 10)
  
#=========================================================================================
# detailed subcluster annotation revealed a few remaining doublets
# individual major cell types were investigated individually in "IVIg2_snRNA_1.1 subcluster annotation.R"
#=========================================================================================
  #load
  Samples.combined<-qread("IVIG_snRNA_harmony_f3.qs")
  
  #get cells removed in subclustering
  cells_to_rm_CM <- readRDS("CM/before f4/cells_to_rm_CM_f4.RDS")
  cells_to_rm_FIB <- readRDS("FIB/before filtering f4/cells_to_rm for FIB f4.RDS")
  cells_to_rm_myl <- readRDS("Myeloid/before filtering f4/cells_to_rm for myeloid f4.RDS")
  cells_to_rm_lym <- readRDS("Lymphoid/before filtering f4/cells_to_rm for lymphoid f4.RDS")
  cells_to_rm_Neuro <- readRDS("NEURO/before filtering f4/cells_to_rm for NEURO f4.RDS")
  cells_to_rm<-c(cells_to_rm_CM,cells_to_rm_FIB,cells_to_rm_myl,cells_to_rm_lym,cells_to_rm_Neuro)

  #remove
  Samples.combined<-subset(Samples.combined,cells=cells_to_rm,invert=T)

  #add subclustering annotation to the integrated data
  Samples.combined.EC<-qread("EC/subclustering of EC_annotated.qs")
  Samples.combined.ADI<-qread("ADI/subclustering of ADI_annotated.qs")
  Samples.combined.NEURO<-qread("NEURO/subclustering of NEURO_annotated.qs")
  Samples.combined.FIB<-qread("FIB/subclustering of FB_annotated.qs")
  Samples.combined.CM<-qread("CM/subclustering of CM_annotated.qs")
  Samples.combined.Lymph<-qread("Lymphoid/subclustering of Lymph_annotated.qs")
  Samples.combined.Myel<-qread("Myeloid/before_f5/subclustering of Myel_annotated.qs")
  Samples.combined.Mural<-qread("Mural/subclustering of Mural_annotated.qs")
  subcluster_meta<-FetchData(Samples.combined.CM,vars = "subcelltypes.short")
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.EC,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.FIB,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.Mural,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.Myel,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.Lymph,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.NEURO,vars = "subcelltypes.short"))
  subcluster_meta<-rbind(subcluster_meta,FetchData(Samples.combined.ADI,vars = "subcelltypes.short"))
  Samples.combined<-AddMetaData(Samples.combined,subcluster_meta)
  Samples.combined$subcelltypes.short<-factor(Samples.combined$subcelltypes.short,levels =   levels(Samples.combined$subcelltypes.short)[c(1:9,12:48,10,11,49:52)])

  DimPlot(Samples.combined,group.by = "subcelltypes.short")

  #save after filter
  qsave(Samples.combined,file = "IVIG_snRNA_harmony_f4.qs")

#=========================================================================================
# cell type specific in depth analysis reveal small amounts of remaining artifacts in the  Myeloid cluster
#=========================================================================================
  Samples.combined.Myel<-qread("subclustering of Myel_annotated_f5.qs")
  #get cell names of Myeloid from overall
  names_Myel <- colnames(Samples.combined)[Samples.combined$celltypes.short %in% c("Myeloid")]
  remaing_contaminations<-names_Myel[!names_Myel %in% c(colnames(Samples.combined.Myel))]
  #rm
  Samples.combined<-subset(Samples.combined,cells=remaing_contaminations,invert=T)
  #adjust order for myeloid
  Samples.combined$subcelltypes.short<-factor(Samples.combined$subcelltypes.short)
  label_order<-1:length(levels(Samples.combined$subcelltypes.short))
  names(label_order)<-levels(Samples.combined$subcelltypes.short)
  label_order[names(label_order) %in% levels(Samples.combined.Myel$subcelltypes.short)]<- label_order[levels(Samples.combined.Myel$subcelltypes.short)]
  label_order<-label_order[label_order]
  Samples.combined$subcelltypes.short<-factor(Samples.combined$subcelltypes.short,levels = names(label_order))

#=========================================================================================
# add meta data category of major celltype annotation but extended to more major cts
#=========================================================================================
  Samples.combined$celltypes.short.ext<-as.character(Samples.combined$celltypes.short)
  #Mural
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("PC_1","PC_2","PC_MT")]<-"PC"
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("SMC_1.1","SMC_1.2","SMC_1.3","SMC_2","proliferating_Mural")]<-"SMC"
  #Myeloid
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("MO_CD16","MO_VCAN")]<-"Monocyte"
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("MAST")]<-"MAST"
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("MP_LYVE1hi_MHCIIlow","MP_LYVE1hi_MHCIIint","MP_LYVE1low_MHCIIhi",
                                                                "MP_FOLR2","MP_TREM2" ,"MP_ISG","MP_proliferating1","MP_proliferating2","MP_NFKB")]<-"Macrophage"
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("cDC")]<-"cDC"
  #Lymphoid
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c("CD4T_naive","CD4T_act","CD4T_reg","CD8T_trans","CD8T_cytox_TEMRA","CD8T_TEM","MAIT_like", "proliferating_Lympho" )]<-"Tcells"
  Samples.combined$celltypes.short.ext[Samples.combined$subcelltypes.short %in% c(  "NK_CD16hi" ,"NK_CD56hi" )]<-"NK"
  #add lvl
  Samples.combined$celltypes.short.ext<-factor(Samples.combined$celltypes.short.ext,levels= c("CM","EC","FIB","PC","SMC","Monocyte","Macrophage","cDC","MAST","Tcells","NK","NEURO","LYMEC","EndoEC","ADI"))
  
  #save
  qsave(Samples.combined,file = "IVIG_snRNA_harmony_f5.qs")
