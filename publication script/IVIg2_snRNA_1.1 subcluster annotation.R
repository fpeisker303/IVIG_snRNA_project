## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#subcluster each major cell type, check filtering and add subtype annotation

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#=========================================================================================
# load f3 filtered snRNA data (doublet, low quality and bad samples already removed)
#=========================================================================================
#load data
sc_IVIg<-qread("IVIG_snRNA_harmony_f3.qs")
Idents(sc_IVIg)<-"celltypes.short"

#=========================================================================================
# load reichart to compare annotation / label transfer
#=========================================================================================
reichart.sc<-qread("public_snRNA_data/Reichart et al snRNA DCM Sience/reichart.sc_meta.data_optimized.qs")
Idents(reichart.sc)<-"cell_type"

#=========================================================================================
# CM
#=========================================================================================
setwd("subclustering/CM/")

#get CM from f3
sc_IVIg.CM <- subset(sc_IVIg,ident="CM")

# rerun harmony in cell type and create new UMAP ---
  sc_IVIg.CM <- sc_IVIg.CM %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
  ElbowPlot(sc_IVIg.CM,reduction = "harmony",ndims = 50)
  
  sc_IVIg.CM <- sc_IVIg.CM %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)

 # transfer label
  reichart.CM<-subset(reichart.sc,ident="cardiac muscle cell")
  sc_IVIg.CM <- RunKNNPredict(
    srt_query = sc_IVIg.CM, srt_ref = reichart.CM,
    ref_group = "cell_states",
    return_full_distance_matrix = TRUE,
    prefix="subcelltype"
  )
  rm(reichart.CM)
  gc()
  
#search for a resonable cluster resolution
  dir.create("cluster_res_opti")
  #cluster resolution loop
  for (res in seq(0.1,1,0.1)) {
    #create folder
    setwd("subclustering/CM/cluster_res_opti")
    dir.create(gsub("[.]","_",as.character(res)))
    setwd(paste0("subclustering/CM/cluster_res_opti/",gsub("[.]","_",as.character(res))))
    #cluster with res
    print(res)
    sc_IVIg.CM <- FindClusters(sc_IVIg.CM, resolution = res)
    #plot UMAP
    DimPlot(sc_IVIg.CM,label = T,label.size = 8,pt.size = 1)+NoLegend()
    ggsave(filename = paste0(res," CM cluster.jpeg"), width=10 , height = 10)

    #plot nfeature and mt pct
    VlnPlot(sc_IVIg.CM,features = "nFeature_RNA",pt.size = 0)+NoLegend()
    ggsave(filename = paste0(res," CM nFeature_RNA vln.jpeg"), width=6 , height = 5)
    VlnPlot(sc_IVIg.CM,features = "percent.mt",pt.size = 0)+NoLegend()
    ggsave(filename = paste0(res," CM percent.mt vln.jpeg"), width=6 , height = 5)
    
    #get marker genes based on pseudobulk DGEA
    sc_IVIg.CM_ps<-create_pseudobulk_object(sc_IVIg.CM,"orig.ident",cell_type_lable = "seurat_clusters",pseudobulk_per_ct = T,  average_meta_feature = c("percent.mt","percent.rb"),cell_number_cutoff = 10)
    sc_IVIg.CM_ps<-pseudobulk.QC.filter(sc_IVIg.CM_ps,min.percent.sample = 0.3) 
    combined_plot<-pseudobulk.QC.plots(sc_IVIg.CM_ps)
    sc_IVIg.CM_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg.CM_ps@assays$RNA@counts)
    sc_IVIg.CM_ps<-NormalizeData_DESeq2(sc_IVIg.CM_ps)
    sc_IVIg.CM_ps<-ScaleData_DESeq2(sc_IVIg.CM_ps)
    
    all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.CM_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = T,filter_pct_exp_per_cluster_cutoff = 0.5)
    all.markers_f<-all.markers[all.markers$logFC>0 & all.markers$adj.P.Val < 0.01,]
    
    top10 <- all.markers[!grepl("^(AC|AP|AL|LINC)\\d+", rownames(all.markers)),] %>% group_by(cluster) %>% top_n(n = 10,wt=t) %>% pull(gene)
    seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.CM$seurat_clusters))))
    
    FeatureHeatmap(sc_IVIg.CM_ps,group.by = "seurat_clusters",features = all.markers_f$gene,feature_split = all.markers_f$cluster,slot="scale.data",
                   features_label = top10,label_size=7,height = 10,width = 2,assay = "DESeq2",
                   group_palcolor = seurat_cluster_colors,feature_split_palcolor =seurat_cluster_colors,
                   anno_terms = T,db = c("KEGG","GO_BP"),topTerm = 10)
    ggsave(filename = paste0(res," CM subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
    
    all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.CM_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = F)#for down regulated no filter by coverage
    all.markers_f_down<-all.markers[all.markers$logFC<0 & all.markers$adj.P.Val < 0.01,]
    top10 <- all.markers[!grepl("^(AC|AP|AL|LINC)\\d+", rownames(all.markers)),] %>% group_by(cluster) %>% top_n(n = 10,wt=-t) %>% pull(gene)
    
    FeatureHeatmap(sc_IVIg.CM_ps,group.by = "seurat_clusters",features = all.markers_f_down$gene,feature_split = all.markers_f_down$cluster,slot="scale.data",
                   features_label = top10,label_size=7,height = 10,width = 2,assay = "DESeq2",
                   group_palcolor = seurat_cluster_colors,feature_split_palcolor =seurat_cluster_colors,
                   anno_terms = F,db = c("KEGG","GO_BP"),topTerm = 10)
    ggsave(filename = paste0(res," CM subcluster allmarkers heatmap grouped_down.jpeg"), width=15 , height = 15,limitsize = F)
    
    #mlm of Reichart marker genes in res (marker gene t-matrix based)
    #make t value matrix
    all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.CM_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = F)
    all.markers_full<-all.markers[all.markers$gene %in% all.markers_f$gene,] #important! needs all.markers without filtering!
    all.markers_full<-dcast(all.markers_full,formula = gene ~ cluster, value.var = "t")
    row.names(all.markers_full)<-all.markers_full$gene
    all.markers_full$gene<-NULL
    all.markers_full<-all.markers_full[!grepl("^(AC|AP|AL)\\d+", rownames(all.markers_full)),]
    
    #Reichart CM marker genes
    reichart_subtype_markers<-as.data.frame(read_xlsx("public_snRNA_data/Reichart et al snRNA DCM Sience/tables of the paper/S6_CM_Marker_Genes_Cell_States.xlsx"))
    sub_states<- unique(sub("_.*", "", colnames(reichart_subtype_markers)))[-1]
    for (type in sub_states) {
      vec<-paste0(type,c("_gene_symbol","_log2_fold-change"))
      df<-reichart_subtype_markers[,vec]
      colnames(df)<-c("gene_symbol","log2_fold-change")
      df$cs<-paste0("Reichart_",type)
      if (type==sub_states[1]) { df_combined<-df}else{df_combined<-rbind(df_combined,df)}}
    df_combined$FC<-2^df_combined$`log2_fold-change`
    # Run mlm
    sample_acts <- decoupleR::run_mlm(mat = all.markers_full,net = df_combined, .source = 'cs', .target = 'gene_symbol',  .mor = 'FC',  minsize = 1)
    # Transform to wide matrix
    sample_acts_mat <- sample_acts %>%
      tidyr::pivot_wider(id_cols = 'condition',names_from = 'source',values_from = 'score') %>%
      tibble::column_to_rownames('condition') %>% as.matrix()
    # Scale per feature
    sample_acts_mat <- scale(sample_acts_mat)
    my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),seq(0.05,2, length.out = floor(100 / 2)))
    #add pval 
    sample_acts <- sample_acts %>% mutate(sig = case_when(p_value < 0.001 ~ "***",p_value < 0.01 ~ "**",p_value < 0.05 ~ "*", TRUE ~ ""))
    # Transform to wide format for annotations, similar to `sample_acts_mat`
    annotations <- sample_acts %>%
      tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',values_from = 'sig') %>%
      tibble::column_to_rownames('condition') %>% as.matrix()
    #plot
    hm<-pheatmap::pheatmap(mat = sample_acts_mat, color = col.ramp(100),cluster_rows = FALSE,cluster_cols = F,border_color = "white",
                           breaks = my_breaks,cellwidth = 20,cellheight = 20,treeheight_row = 20,treeheight_col = 20,
                           display_numbers = annotations,number_color = "black") # Adjust `number_color` as needed
    dev.off()  
    jpeg(filename = paste0("CM Reichart subtype scores based on t-value marker genes ",res,".jpeg"), width=1500 , height =1500,quality = 100,res = 300)
    print(hm)
    dev.off()  
    
  }
  setwd("subclustering/CM/")
  clustree(sc_IVIg.CM)
  ggsave(filename = "CM subcluster clustree.jpeg", width=10 , height = 10)
  
  #remove all clusterings
  meta.data<-sc_IVIg.CM@meta.data
  str_detect(colnames(meta.data),"RNA_snn_res")
  meta.data<-meta.data[,!str_detect(colnames(meta.data),"RNA_snn_res")]
  sc_IVIg.CM@meta.data<-meta.data
  
  #check for remaining doublets
  sc_IVIg.CM<-ScaleData(sc_IVIg.CM,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  FeaturePlotCombined(sc_IVIg.CM,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
  ggsave(filename = "sc_IVIg.CM marker genes.jpeg", width=20 , height = 15)
  
  #annotation in res 0.3
  res=0.3
  sc_IVIg.CM <- FindClusters(sc_IVIg.CM, resolution = res)
  Idents(sc_IVIg.CM)="seurat_clusters"
  types <- list("0"="CM1","1"="CM2","2"="CM3","3"="CM4","4"="CM5","5"="CM6")
  #Rename identities
  sc_IVIg.CM <- RenameIdents(sc_IVIg.CM, types)
  sc_IVIg.CM$subcelltypes.short <- Idents(sc_IVIg.CM)
  cell.levels = c("CM1","CM2","CM3","CM4","CM5","CM6")
  sc_IVIg.CM$subcelltypes.short = factor(sc_IVIg.CM$subcelltypes.short,levels = cell.levels)
  Idents(sc_IVIg.CM)="subcelltypes.short"
  
  #further downstream analysis revealed CM6 as most likly a small proportion of remaining doublets/artifacts
  #save which CM to remove from the total dataset
  cells_to_rm <-WhichCells(sc_IVIg.CM,idents= "CM6")
  saveRDS(cells_to_rm,"cells_to_rm_CM_f4.RDS")  
  
  sc_IVIg.CM<-subset(sc_IVIg.CM,ident="CM6",invert=T)
  sc_IVIg.CM$subcelltypes.short<-factor(sc_IVIg.CM$subcelltypes.short)
  DimPlot(sc_IVIg.CM,label = T,label.size = 8)+NoLegend()
  ggsave(filename = "CM sub res0.3_CM6 rm .jpeg", width=10 , height = 10)
  
  sc_IVIg.CM <- sc_IVIg.CM %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)
  DimPlot(sc_IVIg.CM,label = T,label.size = 8)+NoLegend()
  
  #save with annotation
  qsave(sc_IVIg.CM,file = "./subclustering of CM_annotated.qs")
  
#=========================================================================================
# EC, EndoEC and LymEC
#=========================================================================================
setwd("subclustering/EC")

#get EC from f3
sc_IVIg.EC <- subset(sc_IVIg,ident=c("EC","EndoEC","LYMEC"))

# rerun harmony in cell type and create new UMAP ---
sc_IVIg.EC <- sc_IVIg.EC %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.EC,reduction = "harmony",ndims = 50)

sc_IVIg.EC <- sc_IVIg.EC %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)

# transfer label
reichart.EC<-subset(reichart.sc,ident="endothelial cell")
sc_IVIg.EC <- RunKNNPredict(
  srt_query = sc_IVIg.EC, srt_ref = reichart.EC,
  ref_group = "cell_states",
  return_full_distance_matrix = TRUE,
  prefix="subcelltype"
)
rm(reichart.EC)
gc()

DimPlot(sc_IVIg.EC,group.by = "subcelltype_classification", cols = safe_c)
ggsave(filename = paste0("sc_IVIg.EC KNNPredict_classification of cell_states.jpeg"), width=8 , height = 7)
DimPlot(sc_IVIg.EC,group.by = "subcelltype_classification",split.by = "subcelltype_classification", cols = safe_c,ncol = 3)
ggsave(filename = paste0("sc_IVIg.EC KNNPredict_classification of cell_states_split.jpeg"), width=8 , height = 7)

#markers for the 3 major EC subtypes
#general capillary
capillary.marker<-unique(c(c("CD36","BTNL9","MGLL","ABLIM3"),
                           c("BTNL9","RGCC","CPAMB8","LNX1"),
                           c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","MGLL")))
#artery
artery.marker<-unique(c(c("HEY1","DDK2","NEBL1","EFNB2"),
                        c("NKAIN2","DKK2","TOX","GJA5"),
                        c("SEMA3G","PCSK5","ARL15","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2")))
#vein
vein.marker<-unique(c(c("NRF2F2","ACKR1","TSHZ2","IGFBP5","EPHB4"),
                      c("ACKR1","POSTN","FAM155A","TSHZ2","SLCC2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2"),
                      c("KCNIP4","TPO","TSHZ2","ABCB1")))#KCNIP4-EC
#endocard 1
endo1.marker <- unique(c(c("NRG3","PCDH15","CDH11"),
                         c("TMEM108","INHBA","WNT9B","NRG3")))
#endocard 2
endo2.marker <- unique(c(c("NRG1","PLXNA4","BMP6"),
                         c("NRG1","PPKG1","GATA4","TMEM108")))
#angiogenic
angio.marker<-c("TMEM163","LOXD1","RGCC","NR5A2")

EC.marker<-list(capillary.marker,artery.marker,vein.marker,angio.marker,endo1.marker,endo2.marker)
#score for marker expression
sc_IVIg.EC<-AddModuleScore(sc_IVIg.EC,features = EC.marker,name =c("capillary.marker","artery.marker","vein.marker","angio.marker","endo1.marker","endo2.marker") )

#search for a resonable cluster resolution
dir.create("cluster_res_opti")
#cluster resolution loop
for (res in seq(0.1,0.7,0.1)) {
  #create folder
  setwd("subclustering/EC/cluster_res_opti")
  dir.create(gsub("[.]","_",as.character(res)))
  setwd(paste0("subclustering/EC/cluster_res_opti/",gsub("[.]","_",as.character(res))))
  #cluster with res
  print(res)
  sc_IVIg.EC <- FindClusters(sc_IVIg.EC, resolution = res)
  #plot UMAP
  DimPlot(sc_IVIg.EC,label = T,label.size = 8,pt.size = 1)+NoLegend()
  ggsave(filename = paste0(res," EC cluster.jpeg"), width=10 , height = 10)

  #plot nfeature and mt pct
  VlnPlot(sc_IVIg.EC,features = "nFeature_RNA",pt.size = 0)+NoLegend()
  ggsave(filename = paste0(res," EC nFeature_RNA vln.jpeg"), width=6 , height = 5)
  VlnPlot(sc_IVIg.EC,features = "percent.mt",pt.size = 0)+NoLegend()
  ggsave(filename = paste0(res," EC percent.mt vln.jpeg"), width=6 , height = 5)

  #get marker genes based on pseudobulk DGEA
  sc_IVIg.EC_ps<-create_pseudobulk_object(sc_IVIg.EC,"orig.ident",cell_type_lable = "seurat_clusters",pseudobulk_per_ct = T,  average_meta_feature = c("percent.mt","percent.rb"),cell_number_cutoff = 10)
  sc_IVIg.EC_ps<-pseudobulk.QC.filter(sc_IVIg.EC_ps,min.percent.sample = 0.3)
  combined_plot<-pseudobulk.QC.plots(sc_IVIg.EC_ps)
  sc_IVIg.EC_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg.EC_ps@assays$RNA@counts)
  sc_IVIg.EC_ps<-NormalizeData_DESeq2(sc_IVIg.EC_ps)
  sc_IVIg.EC_ps<-ScaleData_DESeq2(sc_IVIg.EC_ps)

  all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.EC_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = T,filter_pct_exp_per_cluster_cutoff = 0.5)
  all.markers_f<-all.markers[all.markers$logFC>0 & all.markers$adj.P.Val < 0.01,]

  top10 <- all.markers[!grepl("^(AC|AP|AL|LINC)\\d+", rownames(all.markers)),] %>% group_by(cluster) %>% top_n(n = 10,wt=t) %>% pull(gene)
  seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.EC$seurat_clusters))))

  FeatureHeatmap(sc_IVIg.EC_ps,group.by = "seurat_clusters",features = all.markers_f$gene,feature_split = all.markers_f$cluster,slot="scale.data",
                 features_label = top10,label_size=7,height = 10,width = 2,assay = "DESeq2",
                 group_palcolor = seurat_cluster_colors,feature_split_palcolor =seurat_cluster_colors,
                 anno_terms = T,db = c("KEGG","GO_BP"),topTerm = 10)
  ggsave(filename = paste0(res," EC subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)

  all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.EC_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = F)#for down regulated no filter by coverage
  all.markers_f_down<-all.markers[all.markers$logFC<0 & all.markers$adj.P.Val < 0.01,]
  top10 <- all.markers[!grepl("^(AC|AP|AL|LINC)\\d+", rownames(all.markers)),] %>% group_by(cluster) %>% top_n(n = 10,wt=-t) %>% pull(gene)

  FeatureHeatmap(sc_IVIg.EC_ps,group.by = "seurat_clusters",features = all.markers_f_down$gene,feature_split = all.markers_f_down$cluster,slot="scale.data",
                 features_label = top10,label_size=7,height = 10,width = 2,assay = "DESeq2",
                 group_palcolor = seurat_cluster_colors,feature_split_palcolor =seurat_cluster_colors,
                 anno_terms = F,db = c("KEGG","GO_BP"),topTerm = 10)
  ggsave(filename = paste0(res," EC subcluster allmarkers heatmap grouped_down.jpeg"), width=15 , height = 15,limitsize = F)

  #mlm of Reichart marker genes in res (marker gene t-matrix based)
  #make t value matrix
  all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.EC_ps,cluster_identifier = "seurat_clusters",filter_pct_exp_per_cluster = F)
  all.markers_full<-all.markers[all.markers$gene %in% all.markers_f$gene,] #important! needs all.markers without filtering!
  all.markers_full<-dcast(all.markers_full,formula = gene ~ cluster, value.var = "t")
  row.names(all.markers_full)<-all.markers_full$gene
  all.markers_full$gene<-NULL
  all.markers_full<-all.markers_full[!grepl("^(AC|AP|AL)\\d+", rownames(all.markers_full)),]

  #Reichart EC marker genes
  reichart_subtype_markers<-as.data.frame(read_xlsx("public_snRNA_data/Reichart et al snRNA DCM Sience/tables of the paper/S29_EC_Marker_Genes_Cell_States.xlsx"))
  sub_states<- unique(sub("_.*", "", colnames(reichart_subtype_markers)))[-1]
  for (type in sub_states) {
    vec<-paste0(type,c("_gene_symbol","_log2_fold-change"))
    df<-reichart_subtype_markers[,vec]
    colnames(df)<-c("gene_symbol","log2_fold-change")
    df$cs<-paste0("Reichart_",type)
    if (type==sub_states[1]) { df_combined<-df}else{df_combined<-rbind(df_combined,df)}}
  df_combined$FC<-2^df_combined$`log2_fold-change`
  # Run mlm
  sample_acts <- decoupleR::run_mlm(mat = all.markers_full,net = df_combined, .source = 'cs', .target = 'gene_symbol',  .mor = 'FC',  minsize = 1)
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',names_from = 'source',values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>% as.matrix()
  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),seq(0.05,2, length.out = floor(100 / 2)))
  #add pval
  sample_acts <- sample_acts %>% mutate(sig = case_when(p_value < 0.001 ~ "***",p_value < 0.01 ~ "**",p_value < 0.05 ~ "*", TRUE ~ ""))
  # Transform to wide format for annotations, similar to `sample_acts_mat`
  annotations <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',values_from = 'sig') %>%
    tibble::column_to_rownames('condition') %>% as.matrix()
  #plot
  hm<-pheatmap::pheatmap(mat = sample_acts_mat, color = col.ramp(100),cluster_rows = FALSE,cluster_cols = F,border_color = "white",
                         breaks = my_breaks,cellwidth = 20,cellheight = 20,treeheight_row = 20,treeheight_col = 20,
                         display_numbers = annotations,number_color = "black") # Adjust `number_color` as needed
  dev.off()
  jpeg(filename = paste0("EC Reichart subtype scores based on t-value marker genes ",res,".jpeg"), width=1500 , height =1500,quality = 100,res = 300)
  print(hm)
  dev.off()

  #score for marker expression and plot
  VlnPlot(sc_IVIg.EC,features = c("capillary.marker1","artery.marker2","vein.marker3","angio.marker4","endo1.marker5","endo2.marker6"),stack = T,flip = T)
  ggsave(filename = paste0("EC subclustering_EC subtype marker scores.jpeg") , width=20 , height = 30)
}
setwd("subclustering/EC/")
clustree(sc_IVIg.EC)
ggsave(filename = "EC subcluster clustree.jpeg", width=10 , height = 10)

#remove all clusterings
meta.data<-sc_IVIg.EC@meta.data
str_detect(colnames(meta.data),"RNA_snn_res")
meta.data<-meta.data[,!str_detect(colnames(meta.data),"RNA_snn_res")]
sc_IVIg.EC@meta.data<-meta.data

FeaturePlotCombined(sc_IVIg.EC,features = c("capillary.marker1","artery.marker2","vein.marker3","angio.marker4","endo1.marker5","endo2.marker6"))
ggsave(filename = "sc_IVIg.EC endothelial subtype scores.jpeg", width=15 , height = 10)

#check for classic cardiac ct marker genes
sc_IVIg.EC<-ScaleData(sc_IVIg.EC,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
FeaturePlotCombined(sc_IVIg.EC,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
ggsave(filename = "sc_IVIg.EC marker genes.jpeg", width=20 , height = 15)


#annotation in res 0.3
res=0.3
sc_IVIg.EC <- FindClusters(sc_IVIg.EC, resolution = res)
Idents(sc_IVIg.EC)="seurat_clusters"
types <- list("0"="capillaryEC",
              "3"="veinEC",
              "2"="arteryEC",
              "1"="EndoEC",
              "4"="angioEC",
              "5"="LymEC")
#Rename identities
sc_IVIg.EC <- RenameIdents(sc_IVIg.EC, types)
sc_IVIg.EC$subcelltypes.short <- Idents(sc_IVIg.EC)
cell.levels = c("capillaryEC","veinEC","arteryEC","angioEC","EndoEC","LymEC")
sc_IVIg.EC$subcelltypes.short = factor(sc_IVIg.EC$subcelltypes.short,levels = cell.levels)
Idents(sc_IVIg.EC)="subcelltypes.short"

#save with annotation
qsave(sc_IVIg.EC,file = "./subclustering of EC_annotated.qs")


#=========================================================================================
# Fibroblast
#=========================================================================================
setwd("subclustering/FIB/")

#get FIB from f3
sc_IVIg.FIB <- subset(sc_IVIg,ident="FIB")

# rerun harmony in cell type and create new UMAP ---
sc_IVIg.FIB <- sc_IVIg.FIB %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.FIB,reduction = "harmony",ndims = 50)

sc_IVIg.FIB <- sc_IVIg.FIB %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)

# transfer label
reichart.FIB<-subset(reichart.sc,ident="fibroblast of cardiac tissue")
sc_IVIg.FIB <- RunKNNPredict(
  srt_query = sc_IVIg.FIB, srt_ref = reichart.FIB,
  ref_group = "cell_states",
  return_full_distance_matrix = TRUE,
  prefix="subcelltype"
)
rm(reichart.FIB)
gc()

DimPlot(sc_IVIg.FIB,group.by = "subcelltype_classification", cols = safe_c)
ggsave(filename = paste0("sc_IVIg.FIB KNNPredict_classification of cell_states.jpeg"), width=8 , height = 7)
DimPlot(sc_IVIg.FIB,group.by = "subcelltype_classification",split.by = "subcelltype_classification", cols = safe_c,ncol = 3)
ggsave(filename = paste0("sc_IVIg.FIB KNNPredict_classification of cell_states_split.jpeg"), width=8 , height = 7)

#  clustering based annotation
#search for a resonable cluster resolution
plot.list<-NULL
dir.create("cluster_res_heatmaps")
#cluster resolution loop
for (res in seq(0.1,1,0.1)) {
  sc_IVIg.FIB <- FindClusters(sc_IVIg.FIB, resolution = res)
  print(res)
  plot.list[[as.character(res)]]<-DimPlot(sc_IVIg.FIB,label = T,label.size = 8)+NoLegend()

  #get marker genes quick with wilcox test and filter
  all.markers<-FindAllMarkers(sc_IVIg.FIB, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
  all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
  sc_IVIg.FIB<-ScaleData(sc_IVIg.FIB,features = unique(all.markers$gene))
  seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.FIB$seurat_clusters))))
  ht<-GroupHeatmap(sc_IVIg.FIB,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",features_label = top10,label_size=7,height = 10,
                   group_palcolor = seurat_cluster_colors,feature_split_palcolor = seurat_cluster_colors)
  ht
  ggsave(filename = paste0("cluster_res_heatmaps/",res," FIB subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
}
plot_grid(plotlist = plot.list)
ggsave(filename = "FIB subcluster umap seurat cluster resolutions.jpeg", width=20 , height = 15)
clustree(sc_IVIg.FIB)
ggsave(filename = "FIB subcluster clustree.jpeg", width=10 , height = 10)

#classic markers of Fib subtypes
fib.marker <- c("ACTA2","THBS4","SCARA5","POSTN","FN1","COL1A1","COL15A1")
FeaturePlotCombined(sc_IVIg.FIB,features = fib.marker,pt.size = 3)
ggsave(filename = paste0("FIB subclustering_fib.marker.jpeg") , width=30 , height = 30)
DotPlot(sc_IVIg.FIB,features = fib.marker)+coord_flip() + scale_color_gradientn(colours = gradient.col)
ggsave(filename = paste0("FIB subclustering_fib.marker_dp.jpeg") , width=6 , height = 4)

#check for remaining doublets
sc_IVIg.FIB<-ScaleData(sc_IVIg.FIB,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
FeaturePlotCombined(sc_IVIg.FIB,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
ggsave(filename = "sc_IVIg.FIB marker genes.jpeg", width=20 , height = 15)
FeaturePlotCombined(sc_IVIg.FIB,features = c("percent.mt","percent.rb"))
ggsave(filename = "sc_IVIg.FIB percent mt and rb.jpeg", width=10 , height = 5)

#annotation in res 0.5
res=0.5
sc_IVIg.FIB <- FindClusters(sc_IVIg.FIB, resolution = res)
Idents(sc_IVIg.FIB)="seurat_clusters"

DotPlot(sc_IVIg.FIB,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))+coord_flip() + scale_color_gradientn(colours = gradient.col)
ggsave(filename = paste0("sc_IVIg.FIB cardiac ct marker genes.jpeg") , width=6 , height = 4)

#there still seem to be contamination
#remove those remaining
cells_to_rm<-WhichCells(sc_IVIg.FIB,ident=c(6,8))
saveRDS(cells_to_rm,"cells_to_rm for FIB f4.RDS")
sc_IVIg.FIB<-subset(sc_IVIg.FIB,ident=c(6,8),invert=T)

#rerun harmony
sc_IVIg.FIB <- sc_IVIg.FIB %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.FIB,reduction = "harmony",ndims = 50)
sc_IVIg.FIB <- sc_IVIg.FIB %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)
sc_IVIg.FIB <- FindClusters(sc_IVIg.FIB, resolution = res)
Idents(sc_IVIg.FIB)="seurat_clusters"

#classic markers of Fib subtypes
fib.marker <- c("ACTA2","THBS4","SCARA5","POSTN","FN1","COL1A1","COL15A1")
FeaturePlotCombined(sc_IVIg.FIB,features = fib.marker,pt.size = 3)
ggsave(filename = paste0("FIB subclustering_fib.marker_recluster.jpeg") , width=30 , height = 30)
DotPlot(sc_IVIg.FIB,features = fib.marker)+coord_flip() + scale_color_gradientn(colours = gradient.col)
ggsave(filename = paste0("FIB subclustering_fib.marker_dp_recluster.jpeg") , width=6 , height = 4)

#annotation
types <- list("0"="FB1","1"="FB1",
              "2"="FB2",
              "3"="SCARA5_FB",
              "4"="stressed_FB",
              "5"="ECM_FB",
              "6"="Myo_FB")

#Rename identities
sc_IVIg.FIB <- RenameIdents(sc_IVIg.FIB, types)
sc_IVIg.FIB$subcelltypes.short <- Idents(sc_IVIg.FIB)
cell.levels = c("FB1","FB2","stressed_FB","SCARA5_FB","ECM_FB","Myo_FB")
sc_IVIg.FIB$subcelltypes.short = factor(sc_IVIg.FIB$subcelltypes.short,levels = cell.levels)
Idents(sc_IVIg.FIB)="subcelltypes.short"

#save with annotation
qsave(sc_IVIg.FIB,file = "./subclustering of FB_annotated.qs")

#=========================================================================================
# Myeloid immune cells
#=========================================================================================
setwd("subclustering/Myeloid/")

#get Myel from f3
  sc_IVIg.Myel<- subset(sc_IVIg,ident="Myeloid")

# rerun harmony in cell type and create new UMAP ---
  sc_IVIg.Myel <- sc_IVIg.Myel %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony("orig.ident", plot_convergence = TRUE)
  ElbowPlot(sc_IVIg.Myel,reduction = "harmony",ndims = 50)

  sc_IVIg.Myel <- sc_IVIg.Myel %>%
    RunUMAP(reduction = "harmony", dims = 1:25) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25)

# transfer label
  reichart.sc<-qread("public_snRNA_data/Reichart et al snRNA DCM Sience/reichart.sc_meta.data_optimized.qs")
  Idents(reichart.sc)<-"cell_type"
  reichart.Myel<-subset(reichart.sc,ident="myeloid cell")
  sc_IVIg.Myel <- RunKNNPredict(
    srt_query = sc_IVIg.Myel, srt_ref = reichart.Myel,
    ref_group = "cell_states",
    return_full_distance_matrix = TRUE,
    prefix="subcelltype"
  )
  DimPlot(sc_IVIg.Myel,group.by = "subcelltype_classification",cols = dittoColors())
  ggsave(filename = paste0("sc_IVIg.Myel KNNPredict_classification of cell_states.jpeg"), width=8 , height = 7)
  DimPlot(sc_IVIg.Myel,group.by = "subcelltype_classification",split.by = "subcelltype_classification",ncol = 3,cols = dittoColors())
  ggsave(filename = paste0("sc_IVIg.Myel KNNPredict_classification of cell_states_split.jpeg"), width=15 , height = 15)

  #  clustering based annotation
  #search for a resonable cluster resolution
  plot.list<-NULL
  dir.create("cluster_res_heatmaps")
  #cluster resolution loop
  for (res in seq(0.1,1,0.1)) {
    sc_IVIg.Myel <- FindClusters(sc_IVIg.Myel, resolution = res)
    print(res)
    plot.list[[as.character(res)]]<-DimPlot(sc_IVIg.Myel,label = T,label.size = 8)+NoLegend()

    #get marker genes quick with wilcox test and filter
    all.markers<-FindAllMarkers(sc_IVIg.Myel, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
    all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
    top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
    sc_IVIg.Myel<-ScaleData(sc_IVIg.Myel,features = unique(all.markers$gene))
    seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.Myel$seurat_clusters))))
    ht<-GroupHeatmap(sc_IVIg.Myel,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",features_label = top10,label_size=7,height = 10,
                     group_palcolor = seurat_cluster_colors,feature_split_palcolor = seurat_cluster_colors)
    ht
    ggsave(filename = paste0("cluster_res_heatmaps/",res," Myel subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
  }
  plot_grid(plotlist = plot.list)
  ggsave(filename = "Myel subcluster umap seurat cluster resolutions.jpeg", width=20 , height = 15)
  clustree(sc_IVIg.Myel)
  ggsave(filename = "Myel subcluster clustree.jpeg", width=10 , height = 10)

  #annotation in res 0.4
  res=0.4
  sc_IVIg.Myel <- FindClusters(sc_IVIg.Myel, resolution = res)
  Idents(sc_IVIg.Myel)="seurat_clusters"

  #explore gene expression of Myeloid clusters
  dittoBarPlot(sc_IVIg.Myel,var = "subcelltype_classification",group.by = "seurat_clusters",retain.factor.levels = T)
  ggsave(filename = paste0("Myel res 0.4_KNNPredict_classification of cell_states.jpeg") , width=10 , height = 4)
  DotPlot(sc_IVIg.Myel,features = c("LYVE1","HLA-DRA","HLA-DQA1" ,"HLA-DQB1" ,"HLA-DPA1", "HLA-DPB1"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  ggsave(filename = paste0("LYVE1 and MHCII genes.jpeg") , width=7 , height = 4)
  #marker genes of the 3 MP subtypes from Reichart et al supplement
  MP.sub.marker<-c("LYVE1","DAAM2","FGF13","CCDC141","IL2RA","LPAR1","WASHC2A","MAMDC2-AS1","EMP1","SLCO3A1","FRY","CD74","IGSF21","BCL2","HLA-DRB1","USP53")
  DotPlot(subset(sc_IVIg.Myel,ident=c(0,1,2)),features =MP.sub.marker)+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  DotPlot(sc_IVIg.Myel,features = c("TOP2A","MKI67"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#proliferation
  DotPlot(sc_IVIg.Myel,features = c("FCGR3A","VCAN"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#monocyte type
  DotPlot(sc_IVIg.Myel,features = c("MX1","MX2"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#ISG type
  DotPlot(sc_IVIg.Myel,features = c("MYL2","C1QC"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#FOLR2 type
  DotPlot(sc_IVIg.Myel,features = c("TREM2","TPRG1","ITGAX","MYO1E"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#TREM2 type
  DimPlot(sc_IVIg.Myel,label = T,label.size = 8)+NoLegend()
  ggsave(filename = "Myel 0.4 subcluster umap.jpeg", width=8 , height = 8)
  #get res 0.4 marker genes
  all.markers<-FindAllMarkers(sc_IVIg.Myel, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA",test.use = "MAST")
  top5<- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% pull(gene))
  DotPlot(sc_IVIg.Myel, features = rev(top5),dot.scale = 10,assay = "RNA")+ coord_flip()+ scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  ggsave(filename = paste0("Myel res0.4_top5 per cluster dp.jpeg"), width=8 , height = 15)
  #check contaminations
  DotPlot(sc_IVIg.Myel,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","CPA3","GPAM","MMRN1"))+coord_flip() +
      scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(filename = paste0("sc_IVIg.Myel cardiac ct marker genes.jpeg") , width=6 , height = 4)
  DotPlot(sc_IVIg.Myel,features = c("CD3E","CD3D"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#Tcell doublets
  #cluster 10
  DotPlot(sc_IVIg.Myel,features = c("IGHM","CD79B"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#Bcell doublets
  #cluster 12
  DotPlot(sc_IVIg.Myel,features = c("IGKC","IGLC2"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))#Plasma cell doublets
  #cluster 15

  #remove those remaining
  cells_to_rm<-WhichCells(sc_IVIg.Myel,ident=c(10,12,15))
  saveRDS(cells_to_rm,"cells_to_rm for myeloid f4.RDS")
  sc_IVIg.Myel<-subset(sc_IVIg.Myel,ident=c(10,12,15),invert=T)
  sc_IVIg.Myel$seurat_clusters<-factor(sc_IVIg.Myel$seurat_clusters)

  Idents(sc_IVIg.Myel)<-"seurat_clusters"
  #annotation based on Reichart et al 2022 labels
  types <- list("0"="MP_LYVE1hi_MHCIIlow",
                "1"="MP_LYVE1low_MHCIIhi",
                "2"="MP_LYVE1hi_MHCIIint",
                "3"="cDC",
                "4"="MP_FOLR2",
                "5"="MO_CD16",
                "6"="MO_VCAN",
                "7"="MP_TREM2",
                "8"="MAST",
                "9"="MP_ISG",
                "11"="MP_NFKB",
                "13"="MP_proliferating1",
                "14"="MP_proliferating2"
                )
  #Rename identities
  sc_IVIg.Myel <- RenameIdents(sc_IVIg.Myel, types)
  sc_IVIg.Myel$subcelltypes.short <- Idents(sc_IVIg.Myel)
  cell.levels = c("cDC","MO_CD16","MO_VCAN","MP_LYVE1hi_MHCIIlow","MP_LYVE1hi_MHCIIint","MP_LYVE1low_MHCIIhi","MP_FOLR2","MP_TREM2","MP_ISG","MP_NFKB","MP_proliferating1","MP_proliferating2","MAST")
  sc_IVIg.Myel$subcelltypes.short = factor(sc_IVIg.Myel$subcelltypes.short,levels = cell.levels)
  Idents(sc_IVIg.Myel)="subcelltypes.short"
  DimPlot(sc_IVIg.Myel)

  #recalculate the umap
  sc_IVIg.Myel <- sc_IVIg.Myel %>%
    RunUMAP(reduction = "harmony", dims = 1:25) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25)

 #further downstream analysis revealed a remaining cluster with high amounts of CM genes, most likely remaining doublets
  FeaturePlotCombined(sc_IVIg.Myel,features = c("MYH7","MYH6","RYR2","TNN","TNNT2","TNNI3","TPM1"))
  ggsave(filename = paste0("Myel remaing CM doublets.jpeg"), width=20, height = 20)

  cluster_plot_list<-NULL
  for (res in seq(0.1,2,0.1)) {
    sc_IVIg.Myel <- FindClusters(sc_IVIg.Myel, resolution = res)
    cluster_plot_list[[paste0("res_",res)]]<-DimPlot(sc_IVIg.Myel,group.by = "seurat_clusters",label = T)+NoLegend()
  }
  ggarrange(plotlist = cluster_plot_list)
  cluster_plot_list[[5]]
  sc_IVIg.Myel<- FindClusters(sc_IVIg.Myel, resolution = 0.5)
  DimPlot(sc_IVIg.Myel,group.by = "seurat_clusters",label = T)+NoLegend()
  #cluster 8 captures the part of CM doublets
  sc_IVIg.Myel<-subset(sc_IVIg.Myel,ident=8,invert=T)
  DimPlot(sc_IVIg.Myel,group.by = "seurat_clusters",label = T)+NoLegend()
  DimPlot(sc_IVIg.Myel,group.by = "subcelltypes.short",label = T)+NoLegend()

  #recalculate the umap
  sc_IVIg.Myel <- sc_IVIg.Myel %>%
    RunUMAP(reduction = "harmony", dims = 1:25) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25)
  DimPlot(sc_IVIg.Myel,group.by = "subcelltypes.short",label = T)+NoLegend()

  #rm seurat_clustering
  sc_IVIg.Myel@meta.data<-sc_IVIg.Myel@meta.data[,!str_detect(colnames(sc_IVIg.Myel@meta.data),"RNA_snn")]
  Idents(sc_IVIg.Myel)<-"subcelltypes.short"

  #save with annotation
  qsave(sc_IVIg.Myel,file = "./subclustering of Myel_annotated_f5.qs")

#=========================================================================================
# Lymphoid immune cells
#=========================================================================================
  setwd("subclustering/Lymphoid/")

#get Lymph from f3
  sc_IVIg.Lymph<- subset(sc_IVIg,ident="Lymphoid")

# rerun harmony in cell type and create new UMAP ---
  sc_IVIg.Lymph <- sc_IVIg.Lymph %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony("orig.ident", plot_convergence = TRUE)
  ElbowPlot(sc_IVIg.Lymph,reduction = "harmony",ndims = 50)

  sc_IVIg.Lymph <- sc_IVIg.Lymph %>%
    RunUMAP(reduction = "harmony", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20)

# transfer label
    #using reichart
    reichart.sc<-qread("public_snRNA_data/Reichart et al snRNA DCM Sience/reichart.sc_meta.data_optimized.qs")
    Idents(reichart.sc)<-"cell_type"
    reichart.Lymph<-subset(reichart.sc,ident="lymphocyte")
    sc_IVIg.Lymph <- RunKNNPredict(
      srt_query = sc_IVIg.Lymph, srt_ref = reichart.Lymph,
      ref_group = "cell_states",
      return_full_distance_matrix = TRUE,
      prefix="subcelltype"
    )
    DimPlot(sc_IVIg.Lymph,group.by = "subcelltype_classification",cols = dittoColors())
    ggsave(filename = paste0("sc_IVIg.Lymph KNNPredict_classification of cell_states.jpeg"), width=8 , height = 7)
    DimPlot(sc_IVIg.Lymph,group.by = "subcelltype_classification",split.by = "subcelltype_classification",ncol = 3,cols = dittoColors())
    ggsave(filename = paste0("sc_IVIg.Lymph KNNPredict_classification of cell_states_split.jpeg"), width=15 , height = 15)

    #using Zhu2022
    Zhu2022.sc<-qread("public_snRNA_data/Zhu 2022 PBMC after checkpoint inhib_myocarditis/Zhu2022_snRNA_harmony_all_filtered.qs")
    Idents(Zhu2022.sc)<-"celltypes.short_l1"
    Zhu2022.Lymph<-subset(Zhu2022.sc,ident=c("Tcells_CD4", "Tcells_CD8", "NK", "Bcells"))
    DefaultAssay(Zhu2022.Lymph)<-"RNA"
    Zhu2022.Lymph<-NormalizeData(Zhu2022.Lymph)
    sc_IVIg.Lymph <- RunKNNPredict(
      srt_query = sc_IVIg.Lymph, srt_ref = Zhu2022.Lymph,
      ref_group = "celltypes.short_l2",
      return_full_distance_matrix = TRUE,
      prefix="Zhu_subcelltype"
    )
    DimPlot(sc_IVIg.Lymph,group.by = "Zhu_subcelltype_classification",cols = dittoColors())
    ggsave(filename = paste0("sc_IVIg.Lymph KNNPredict_Zhu_subcelltype_classification of cell_states.jpeg"), width=8 , height = 7)

  #  clustering based annotation
  #search for a resonable cluster resolution
  plot.list<-NULL
  dir.create("cluster_res_heatmaps")
  #cluster resolution loop
  for (res in seq(0.1,1,0.1)) {
    sc_IVIg.Lymph <- FindClusters(sc_IVIg.Lymph, resolution = res)
    print(res)
    plot.list[[as.character(res)]]<-DimPlot(sc_IVIg.Lymph,label = T,label.size = 8)+NoLegend()
 
    #get marker genes quick with wilcox test and filter
    all.markers<-FindAllMarkers(sc_IVIg.Lymph, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
    all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
    top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
    sc_IVIg.Lymph<-ScaleData(sc_IVIg.Lymph,features = unique(all.markers$gene))
    seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.Lymph$seurat_clusters))))
    ht<-GroupHeatmap(sc_IVIg.Lymph,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",features_label = top10,label_size=7,height = 10,
                     group_palcolor = seurat_cluster_colors,feature_split_palcolor = seurat_cluster_colors)
    ht
    ggsave(filename = paste0("cluster_res_heatmaps/",res," Lymph subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
  }
  plot_grid(plotlist = plot.list)
  ggsave(filename = "Lymph subcluster umap seurat cluster resolutions.jpeg", width=20 , height = 15)
  clustree(sc_IVIg.Lymph)
  ggsave(filename = "Lymph subcluster clustree.jpeg", width=10 , height = 10)

  #check mt% and rb%
  FeaturePlotCombined(sc_IVIg.Lymph,features = c("percent.mt","percent.rb"))
  ggsave(filename = "sc_IVIg.Lymph percent mt and rb.jpeg", width=10 , height = 5)
  FeaturePlotCombined(sc_IVIg.Lymph,features = c("nFeature_RNA","nCount_RNA"))

  #check contaminations
  DotPlot(sc_IVIg.Lymph,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","CPA3","GPAM","MMRN1"))+coord_flip() +
    scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(filename = paste0("sc_IVIg.Lymph cardiac ct marker genes.jpeg") , width=6 , height =4)

  #annotation in res 0.8
  res=0.8
  sc_IVIg.Lymph <- FindClusters(sc_IVIg.Lymph, resolution = res)
  Idents(sc_IVIg.Lymph)="seurat_clusters"

  #explore gene expression of Lymphoid clusters
  dittoBarPlot(sc_IVIg.Lymph,var = "subcelltype_classification",group.by = "seurat_clusters",retain.factor.levels = T)
  ggsave(filename = "Lymph res0.8 predicted KNN.jpeg", width=8 , height = 8)
  DotPlot(sc_IVIg.Lymph,features = c("CD3E","CD4","LEF1","CD69","FOXP3","CD8A","GZMB","GZMK","TIGIT","KLRB1","KLRF1","NCAM1","CCR6","MKI67","TBX21","FCGR3A","KIT"))+
    coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  DotPlot(sc_IVIg.Lymph,features = c("TOP2A","MKI67"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  DimPlot(sc_IVIg.Lymph,label = T,label.size = 8)+NoLegend()
  ggsave(filename = "Lymph res0.8 subcluster umap.jpeg", width=8 , height = 8)

  #Tcell subtype marker genes
  #Th subset marker (based on perplexcity, maybe not optimal ref!)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("TBX21","IFNG") ,reduction = "umap", ncol = 3)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("GATA3","IL4","IL5","IL10","IL13") ,reduction = "umap", ncol = 3)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("RORC","IL17","IL22") ,reduction = "umap", ncol = 3)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("BCL6","IL21") ,reduction = "umap", ncol = 3)
  #TRM
  FeatureDimPlot(sc_IVIg.Lymph, features = c("CXCR6","ITGA1","IFNG","ITGAE","CCL4","CCL3") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features = c("CD4","CD8A","CXCR6","ITGA1","IFNG","ITGAE","CCL4","CCL3"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #CD4 TCM (activation)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("CCR7","SELL","TCF7","IL2","TNF","IL4R") ,reduction = "umap", ncol = 3)
  DotPlot(sc_IVIg.Lymph,features = c("CD4","CD8A","CCR7","SELL","TCF7","IL2","TNF","IL4R"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #Treg
  FeatureDimPlot(sc_IVIg.Lymph, features = c("FOXP3","IL2RA","CTLA4") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features = c("CD4","CD8A","FOXP3","IL2RA","CTLA4"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #CD8 TEM
  FeatureDimPlot(sc_IVIg.Lymph, features = c("CCL5","GZMK","GZMB") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features = c("CD4","CD8A","CCL5","GZMK","GZMB") )+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #CD8 TEM (tissue)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("TUBA1A","TUBB","VIM","LGALS1") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features =  c("CD4","CD8A","TUBA1A","TUBB","VIM","LGALS1"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #CD8 TEM (activated)
  FeatureDimPlot(sc_IVIg.Lymph, features = c("CCL5","GZMK","GZMB") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features =c("CD4","CD8A","CCL5","GZMK","GZMB") )+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  #CD8 TEMRA
  FeatureDimPlot(sc_IVIg.Lymph, features = c("PTPRC","PRF1","NKG7","HOPX","IFIT3","KLRD1","MYO1F","CCL5","GZMH") ,reduction = "umap", ncol = 3)
  DotPlot(subset(sc_IVIg.Lymph,ident=c(1,2,3,5,6,9,10)),features = c("CD4","CD8A","PTPRC","PRF1","NKG7","HOPX","IFIT3","KLRD1","MYO1F","CCL5","GZMH"))+coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))

  #cluster 4 and 8 are very unclear
  sc_IVIg.Lymph_sub<-subset(sc_IVIg.Lymph,ident=c(4,8))
  sc_IVIg.Lymph_sub <- FindClusters(sc_IVIg.Lymph_sub, resolution = 1)
  DimPlot(sc_IVIg.Lymph_sub,label = T,label.size = 8)+NoLegend()
  ggsave(filename = "Lymph cl4 and cl8 sub res1.jpeg", width=8 , height = 8)
  DotPlot(sc_IVIg.Lymph_sub,features = c("CD3E","CD4","LEF1","CD69","FOXP3","CD8A","GZMB","CD16","KLRB1","KLRF1","NCAM1","FCGR3A","CCR6","MKI67","TBX21","KIT"))+
    coord_flip() + scale_color_gradientn(colours = gradient.col)+theme(axis.text.x = element_text(angle = 45,hjust =1))
  all.markers<-FindAllMarkers(sc_IVIg.Lymph_sub, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA",test.use = "MAST")
  #result of cl4 and cl8 subcluster check:
  #4 = NK_CD56 high (has no CD16/FCGR3A)
  #all remaining cells are unclear and seem high in %mt. also no abundance difference between IVIG and PLACEBO. there more removed as unclear
  cells_to_rm<-WhichCells(sc_IVIg.Lymph_sub,ident=4,invert=T)

  #also cl7 in the NK look like high %mt
  cells_to_rm<-c(cells_to_rm,WhichCells(sc_IVIg.Lymph,ident=7))
  saveRDS(cells_to_rm,"cells_to_rm for lymphoid f4.RDS")
  #remove those remaning
  sc_IVIg.Lymph<-subset(sc_IVIg.Lymph,cells=cells_to_rm,invert=T)
  sc_IVIg.Lymph$seurat_clusters<-factor(sc_IVIg.Lymph$seurat_clusters)

  Idents(sc_IVIg.Lymph)<-"seurat_clusters"
  #annotation based on Reichart et al 2022 labels
  types <- list("0"="NK_CD16hi",
                "1"="CD4T_act",
                "2"="CD8T_cytox_TEMRA",
                "3"="CD8T_trans",
                "5"="CD8T_TEM",
                "6"="CD4T_naive",
                "8"="NK_CD56hi",
                "9"="CD4T_reg",
                "10"="MAIT_like",
                "11"="proliferating_Lympho"
                )
  #Rename identities
  sc_IVIg.Lymph <- RenameIdents(sc_IVIg.Lymph, types)
  sc_IVIg.Lymph$subcelltypes.short <- Idents(sc_IVIg.Lymph)
  cell.levels = c("CD4T_naive","CD4T_act","CD4T_reg","CD8T_trans","CD8T_cytox_TEMRA","CD8T_TEM","MAIT_like","NK_CD16hi","NK_CD56hi","proliferating_Lympho")
  sc_IVIg.Lymph$subcelltypes.short = factor(sc_IVIg.Lymph$subcelltypes.short,levels = cell.levels)
  Idents(sc_IVIg.Lymph)="subcelltypes.short"
  DimPlot(sc_IVIg.Lymph)

  #save with annotation
  qsave(sc_IVIg.Lymph,file = "./subclustering of Lymph_annotated.qs")

#=========================================================================================
# Mural
#=========================================================================================
setwd("subclustering/Mural/")

#get Mural from f3
sc_IVIg.Mural <- subset(sc_IVIg,ident="Mural")

# rerun harmony in cell type and create new UMAP ---
sc_IVIg.Mural <- sc_IVIg.Mural %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.Mural,reduction = "harmony",ndims = 50)

sc_IVIg.Mural <- sc_IVIg.Mural %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)


#  clustering based annotation
#search for a resonable cluster resolution
plot.list<-NULL
dir.create("cluster_res_heatmaps")
#cluster resolution loop
for (res in seq(0.1,1,0.1)) {
  sc_IVIg.Mural <- FindClusters(sc_IVIg.Mural, resolution = res)
  print(res)
  plot.list[[as.character(res)]]<-DimPlot(sc_IVIg.Mural,label = T,label.size = 8)+NoLegend()

  #get marker genes quick with wilcox test and filter
  all.markers<-FindAllMarkers(sc_IVIg.Mural, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
  all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
  sc_IVIg.Mural<-ScaleData(sc_IVIg.Mural,features = unique(all.markers$gene))
  seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.Mural$seurat_clusters))))
  ht<-GroupHeatmap(sc_IVIg.Mural,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",features_label = top10,label_size=7,height = 10,
                   group_palcolor = seurat_cluster_colors,feature_split_palcolor = seurat_cluster_colors)
  ht
  ggsave(filename = paste0("cluster_res_heatmaps/",res," Mural subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
}
plot_grid(plotlist = plot.list)
ggsave(filename = "Mural subcluster umap seurat cluster resolutions.jpeg", width=20 , height = 15)
clustree(sc_IVIg.Mural)
ggsave(filename = "Mural subcluster clustree.jpeg", width=10 , height = 10)

#check mt% and rb%
FeaturePlotCombined(sc_IVIg.Mural,features = c("percent.mt","percent.rb"))
ggsave(filename = "sc_IVIg.Mural percent mt and rb.jpeg", width=10 , height = 5)
FeaturePlotCombined(sc_IVIg.Mural,features = c("nFeature_RNA","nCount_RNA"))

#check contaminations
DotPlot(sc_IVIg.Mural,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","CPA3","GPAM","MMRN1"))+coord_flip() +
  scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0("sc_IVIg.Mural cardiac ct marker genes.jpeg") , width=6 , height = 4)
sc_IVIg.Mural<-ScaleData(sc_IVIg.Mural,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
FeaturePlotCombined(sc_IVIg.Mural,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
ggsave(filename = "sc_IVIg.Mural marker genes.jpeg", width=20 , height = 15)

all.markers<-FindAllMarkers(sc_IVIg.Mural, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA",test.use = "MAST")
all.markers.down<-all.markers[all.markers$avg_log2FC<0 & all.markers$p_val_adj<0.05,]
all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
sc_IVIg.Mural<-ScaleData(sc_IVIg.Mural,features = unique(all.markers$gene))

ht<-GroupHeatmap(sc_IVIg.Mural,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",
                 features_label = top10,label_size=7,height = 10,width = 2,
                 anno_terms = T,anno_features = T,db = c("KEGG","GO_BP"),topTerm = 10)
ht
ggsave(filename = paste0("Mural res0.4 allmarkers heatmap grouped.jpeg"), width=20 , height = 20,limitsize = F)

DimPlot(sc_IVIg.Mural,label = T,label.size = 10)

DotPlot(sc_IVIg.Mural,features = c("MKI67","TOP2A"))+coord_flip() + scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
FeaturePlotCombined(sc_IVIg.Mural,features = c("RGS5","MYH11","CNN1","ACTA2","ABCC9","KCNJ8","IGF1R"),addClustPlot = T)

# transfer label
#using reichart
reichart.sc<-qread("public_snRNA_data/Reichart et al snRNA DCM Sience/reichart.sc_meta.data_optimized.qs")
Idents(reichart.sc)<-"cell_type"
reichart.Mural<-subset(reichart.sc,ident="mural cell")
sc_IVIg.Mural <- RunKNNPredict(
  srt_query = sc_IVIg.Mural, srt_ref = reichart.Mural,
  ref_group = "cell_states",
  return_full_distance_matrix = TRUE,
  prefix="subcelltype"
)
DimPlot(sc_IVIg.Mural,group.by = "subcelltype_classification",cols = dittoColors())
ggsave(filename = paste0("sc_IVIg.Mural KNNPredict_classification of cell_states.jpeg"), width=8 , height = 7)
DimPlot(sc_IVIg.Mural,group.by = "subcelltype_classification",split.by = "subcelltype_classification",ncol = 3,cols = dittoColors())
ggsave(filename = paste0("sc_IVIg.Mural KNNPredict_classification of cell_states_split.jpeg"), width=15 , height = 15)

#annotation in res 0.4
res=0.4
sc_IVIg.Mural <- FindClusters(sc_IVIg.Mural, resolution = res)
Idents(sc_IVIg.Mural)="seurat_clusters"

#annotation based on Reichart et al 2022 labels
types <- list("0"="PC_1",
              "1"="PC_1",
              "2"="SMC_1.1",
              "3"="SMC_1.2",
              "4"="PC_MT",
              "5"="SMC_2",
              "6"="PC_2",
              "7"="SMC_1.3",
              "8"="proliferating_Mural"
)
#Rename identities
sc_IVIg.Mural <- RenameIdents(sc_IVIg.Mural, types)
sc_IVIg.Mural$subcelltypes.short <- Idents(sc_IVIg.Mural)
cell.levels = c("PC_1","PC_2","PC_MT","SMC_1.1","SMC_1.2","SMC_1.3","SMC_2","proliferating_Mural")
sc_IVIg.Mural$subcelltypes.short = factor(sc_IVIg.Mural$subcelltypes.short,levels = cell.levels)
Idents(sc_IVIg.Mural)="subcelltypes.short"
DimPlot(sc_IVIg.Mural)

#save with annotation
qsave(sc_IVIg.Mural,file = "./subclustering of Mural_annotated.qs")

#=========================================================================================
# NEURO
#=========================================================================================
setwd("subclustering/NEURO/")

#get NEURO from f3
sc_IVIg.NEURO <- subset(sc_IVIg,ident="NEURO")

# rerun harmony in cell type and create new UMAP ---
sc_IVIg.NEURO <- sc_IVIg.NEURO %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.NEURO,reduction = "harmony",ndims = 50)

sc_IVIg.NEURO <- sc_IVIg.NEURO %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15)

#  clustering based annotation
#search for a resonable cluster resolution
plot.list<-NULL
dir.create("cluster_res_heatmaps")
#cluster resolution loop
for (res in seq(0.1,0.5,0.1)) {
  sc_IVIg.NEURO <- FindClusters(sc_IVIg.NEURO, resolution = res)
  print(res)
  plot.list[[as.character(res)]]<-DimPlot(sc_IVIg.NEURO,label = T,label.size = 8)+NoLegend()
 
  #get marker genes quick with wilcox test and filter
  all.markers<-FindAllMarkers(sc_IVIg.NEURO, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
  all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
  sc_IVIg.NEURO<-ScaleData(sc_IVIg.NEURO,features = unique(all.markers$gene))
  seurat_cluster_colors<-list(gg_color_hue(length(levels(sc_IVIg.NEURO$seurat_clusters))))
  ht<-GroupHeatmap(sc_IVIg.NEURO,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",features_label = top10,label_size=7,height = 10,
                   group_palcolor = seurat_cluster_colors,feature_split_palcolor = seurat_cluster_colors)
  ht
  ggsave(filename = paste0("cluster_res_heatmaps/",res," NEURO subcluster allmarkers heatmap grouped.jpeg"), width=15 , height = 15,limitsize = F)
}
plot_grid(plotlist = plot.list)
ggsave(filename = "NEURO subcluster umap seurat cluster resolutions.jpeg", width=20 , height = 15)
sc_IVIg.NEURO$RNA_snn_res.1<-NULL
clustree(sc_IVIg.NEURO)
ggsave(filename = "NEURO subcluster clustree.jpeg", width=10 , height = 10)

#check mt% and rb%
FeaturePlotCombined(sc_IVIg.NEURO,features = c("percent.mt","percent.rb"))
ggsave(filename = "sc_IVIg.NEURO percent mt and rb.jpeg", width=10 , height = 5)
FeaturePlotCombined(sc_IVIg.NEURO,features = c("nFeature_RNA","nCount_RNA"))

#check contaminations
DotPlot(sc_IVIg.NEURO,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","CPA3","GPAM","MMRN1"))+coord_flip() +
  scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0("sc_IVIg.NEURO cardiac ct marker genes.jpeg") , width=6 , height = 4)
sc_IVIg.NEURO<-ScaleData(sc_IVIg.NEURO,features =   c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
FeaturePlotCombined(sc_IVIg.NEURO,features = c("RYR2","VWF","DCN","PTPRC","RGS5","MYH11","NRXN1","MKI67","KIT","GPAM","MMRN1"))
ggsave(filename = "sc_IVIg.NEURO marker genes.jpeg", width=20 , height = 15)

all.markers<-FindAllMarkers(sc_IVIg.NEURO, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA",test.use = "MAST")
all.markers.down<-all.markers[all.markers$avg_log2FC<0 & all.markers$p_val_adj<0.05,]
all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
sc_IVIg.NEURO<-ScaleData(sc_IVIg.NEURO,features = unique(all.markers$gene))

ht<-GroupHeatmap(sc_IVIg.NEURO,group.by = "seurat_clusters",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",
                 features_label = top10,label_size=7,height = 10,width = 2,
                 anno_terms = T,anno_features = T,db = c("KEGG","GO_BP"),topTerm = 10)
ht
ggsave(filename = paste0("NEURO res0.4 allmarkers heatmap grouped.jpeg"), width=20 , height = 20,limitsize = F)

DimPlot(sc_IVIg.NEURO,label = T,label.size = 10)

DotPlot(sc_IVIg.NEURO,features = c("MKI67","TOP2A"))+coord_flip() + scale_color_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(angle = 45,hjust = 1))

#two clusters seem to be doublets
#in res 0.1 cl 1 and 2
#remove those remaining
res=0.1
sc_IVIg.NEURO <- FindClusters(sc_IVIg.NEURO, resolution = res)
cells_to_rm<-WhichCells(sc_IVIg.NEURO,ident=c(1,2))
saveRDS(cells_to_rm,"cells_to_rm for NEURO f4.RDS")
sc_IVIg.NEURO<-subset(sc_IVIg.NEURO,ident=c(1,2),invert=T)

#annotation in res 0.4
res=0.2
sc_IVIg.NEURO <- FindClusters(sc_IVIg.NEURO, resolution = res)
Idents(sc_IVIg.NEURO)="seurat_clusters"

#annotation based on Reichart et al 2022 labels
types <- list("0"="NEURO_1",
              "1"="NEURO_2",
              "2"="NEURO_3"
)
#Rename identities
sc_IVIg.NEURO <- RenameIdents(sc_IVIg.NEURO, types)
sc_IVIg.NEURO$subcelltypes.short <- Idents(sc_IVIg.NEURO)
cell.levels = c("NEURO_1","NEURO_2","NEURO_3")
sc_IVIg.NEURO$subcelltypes.short = factor(sc_IVIg.NEURO$subcelltypes.short,levels = cell.levels)
Idents(sc_IVIg.NEURO)="subcelltypes.short"
DimPlot(sc_IVIg.NEURO)

#save with annotation
qsave(sc_IVIg.NEURO,file = "./subclustering of NEURO_annotated.qs")

#=========================================================================================
# ADI
#=========================================================================================
setwd("subclustering/ADI/")

#get ADI from f3
sc_IVIg.ADI <- subset(sc_IVIg,ident="ADI")

# rerun harmony in cell type and create new UMAP ---
sc_IVIg.ADI <- sc_IVIg.ADI %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
ElbowPlot(sc_IVIg.ADI,reduction = "harmony",ndims = 50)

sc_IVIg.ADI <- sc_IVIg.ADI %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15)

sc_IVIg.ADI$subcelltypes.short<-as.character(sc_IVIg.ADI$celltypes.short)

#save with annotation
qsave(sc_IVIg.ADI,file = "./subclustering of ADI_annotated.qs")
