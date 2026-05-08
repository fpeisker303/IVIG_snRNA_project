## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#subcluster analysis for CM cells

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#load
sc_IVIg.CM<-qread("subclustering of CM_annotated.qs")

#======================================================
# plots
#======================================================
  DimPlot(sc_IVIg.CM,label = T,label.size = 8,cols = col.CM)+NoLegend()
  ggsave(filename = "CM subcluster umap annotation.jpeg", width=8 , height = 8)
  
#======================================================
# marker genes
#======================================================
  Idents(sc_IVIg.CM)<-"subcelltypes.short"
  all.markers<-FindAllMarkers(sc_IVIg.CM, min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA",test.use = "MAST")
  # qsave(all.markers,"all.markers.qs")
  # all.markers<-qread("./all.markers.qs")
  all.markers<-all.markers[all.markers$avg_log2FC>0 & all.markers$p_val_adj<0.05,]
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10,wt=avg_log2FC) %>% pull(gene)
  sc_IVIg.CM<-ScaleData(sc_IVIg.CM,features = unique(all.markers$gene))

  ht<-GroupHeatmap(sc_IVIg.CM,group.by = "subcelltypes.short",features = all.markers$gene,feature_split = all.markers$cluster,slot="scale.data",
                   features_label = top10,label_size=7,height = 10,width = 2,
                   group_palcolor = list(col.CM),feature_split_palcolor = list(col.CM))
  ht
  ggsave(filename = paste0("CM subcluster seurat clusters res 0.4 allmarkers heatmap grouped.jpeg"), width=20 , height = 20,limitsize = F)

#GO and KEGG with marker genes
  all.markers<-qread("./all.markers.qs")
  all.markers<-all.markers[!grepl("^(AC|AP|AL)\\d+", all.markers$gene),]
  marker_KEGG<-get.markergenes.enrichr(sc.object = sc_IVIg.CM,sample = "CM",test = "MAST",all.markers = all.markers,clusters = "subcelltypes.short",dbs = "KEGG_2021_Human")
  plot.markergenes.enrichr(sample = "CM",go.results = marker_KEGG,filter.disease.KEGG = T,custom.h = 6.6,custom.w = 5,angle = 45,hjust = 1)
  plot.markergenes.GOs(sample = "CM",go.results = marker_GO,custom.h = 9.5,custom.w = 10,angle = 45,hjust = 1)


#======================================================
# composition analysis
#======================================================
  #boxplot of logFC between PRE and POST, calculated for each patient individually
  composition.ratio.box.plot(sc_IVIg.CM,group1<-"subcelltypes.short", group2<-"condition",sampleIDs = "orig.ident",cols = col.ivig.condi[c(2,4)],add_statistic = F,capY = 3)
  ggsave(filename = "CM subcluster ratio change of sample pair_no stats.jpeg", width=8 , height = 8)
  ggsave(filename = "CM subcluster ratio change of sample pair_no stats.svg", width=8 , height = 8)
  
#======================================================
# mlm on t-values of subcluster markers (Progeny, KEGG, kuppe, kanemaru, Reichart)
#======================================================

#get cluster marker based on pseudobulking
sc_IVIg.CM_ps<-create_pseudobulk_object(sc_IVIg.CM,"orig.ident",cell_type_lable = "subcelltypes.short",pseudobulk_per_ct = T,  average_meta_feature = c("percent.mt","percent.rb"))
sc_IVIg.CM_ps<-pseudobulk.QC.filter(sc_IVIg.CM_ps,min.percent.sample = 0.5)
combined_plot<-pseudobulk.QC.plots(sc_IVIg.CM_ps)
sc_IVIg.CM_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg.CM_ps@assays$RNA@counts)
sc_IVIg.CM_ps<-NormalizeData_DESeq2(sc_IVIg.CM_ps)
sc_IVIg.CM_ps<-ScaleData_DESeq2(sc_IVIg.CM_ps)

all.markers<-FindAllMarkers_pseudobulk(sc_IVIg.CM_ps,cluster_identifier = "subcelltypes_short")
all.markers_f<-all.markers[all.markers$logFC>0 & all.markers$adj.P.Val < 0.01,]

#make t value matrix
all.markers_full<-all.markers[all.markers$gene %in% all.markers_f$gene,]
all.markers_full<-dcast(all.markers_full,formula = gene ~ cluster, value.var = "t")
row.names(all.markers_full)<-all.markers_full$gene
all.markers_full$gene<-NULL
all.markers_full<-all.markers_full[!grepl("^(AC|AP|AL)\\d+", rownames(all.markers_full)),]


#score expression of Kuppe CM marker genes in marker genes
  file_path <- "public_snRNA_data/CKuppe_MI_2022/Supplementary Table 10  Marker genes for cardiomyocytes states.xlsx"
  sheet_names <- excel_sheets(file_path)
  
  #use FC as weight for mlm
  for (sheet in sheet_names) {
    df.sheet<-as.data.frame(read_excel(file_path, sheet = sheet))
    df.sheet<-arrange(df.sheet,-avg_log2FC)
    df.sheet$FC<-2^df.sheet$avg_log2FC
    df.sheet$label<-sheet
    df<-df.sheet[c("gene","FC","label")]
    if (sheet==sheet_names[1]) {
      combined_df<-df
    }else{
      combined_df<-rbind(combined_df,df)
    }
  }
  combined_df$label[combined_df$label=="vCM2"]<-paste0(combined_df$label[combined_df$label=="vCM2"],"_pre_stressed")
  combined_df$label[combined_df$label=="vCM3"]<-paste0(combined_df$label[combined_df$label=="vCM3"],"_stressed")

  # Run mlm
  sample_acts <- decoupleR::run_mlm(mat = all.markers_full,
                                    net = combined_df,
                                    .source = 'label',
                                    .target = 'gene',
                                    .mor = 'FC',
                                    minsize = 1)
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)

  # Color scale
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)

  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05,2, length.out = floor(100 / 2)))
  #add pval
  sample_acts <- sample_acts %>%
    mutate(sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  # Transform to wide format for annotations, similar to `sample_acts_mat`
  annotations <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'sig') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  #plot
  hm<-pheatmap::pheatmap(mat = sample_acts_mat,
                         color = colors.use,
                         cluster_rows = FALSE,
                         cluster_cols = F,
                         border_color = "white",
                         breaks = my_breaks,
                         cellwidth = 20,
                         cellheight = 20,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         display_numbers = annotations,
                         number_color = "black") # Adjust `number_color` as needed
  jpeg(filename = paste0("Kuppe CMsubtype scores based on t-value marker genes_all CM sub.jpeg"), width=1500 , height =1500,quality = 100,res = 300)
  print(hm)
  dev.off()
  svg(filename = paste0("Kuppe CMsubtype scores based on t-value marker genes_all CM sub.svg"), width=4,height = 4)
  print(hm)
  dev.off()

#score expression of Kanemaru CM marker genes in marker genes
  all.markers_kanemaruCM<-qread("public_snRNA_data/Kanemaru_Heart_cell_atlas_2023/all.markers_only Kanemaru CMsubstates.qs")
  all.markers_kanemaruCM<-all.markers_kanemaruCM[all.markers_kanemaruCM$avg_log2FC>0 & all.markers_kanemaruCM$p_val_adj<0.05,]
  table(all.markers_kanemaruCM$cluster)
  all.markers_kanemaruCM$cluster<-paste0("Kanemaru_",all.markers_kanemaruCM$cluster)

  #use lof2FC as weight for mlm
  for (type in unique(all.markers_kanemaruCM$cluster)) {
    if (type==unique(all.markers_kanemaruCM$cluster)[1]) {
      df_combined<-all.markers_kanemaruCM[all.markers_kanemaruCM$cluster==type,c("cluster","avg_log2FC","gene")]
    }else{
      df_combined<-rbind(df_combined,all.markers_kanemaruCM[all.markers_kanemaruCM$cluster==type,c("cluster","avg_log2FC","gene")])
    }
  }
  df_combined$FC<-2^df_combined$avg_log2FC

  # Run mlm
  sample_acts <- decoupleR::run_mlm(mat = all.markers_full,
                                    net = df_combined,
                                    .source = 'cluster',
                                    .target = 'gene',
                                    .mor = 'FC',
                                    minsize = 1)
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)

  # Color scale
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)

  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05,2, length.out = floor(100 / 2)))
  #add pval
  sample_acts <- sample_acts %>%
    mutate(sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  # Transform to wide format for annotations, similar to `sample_acts_mat`
  annotations <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'sig') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  #plot
  hm<-pheatmap::pheatmap(mat = sample_acts_mat,
                         color = colors.use,
                         cluster_rows = FALSE,
                         cluster_cols = F,
                         border_color = "white",
                         breaks = my_breaks,
                         cellwidth = 20,
                         cellheight = 20,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         display_numbers = annotations,
                         number_color = "black") # Adjust `number_color` as needed
  jpeg(filename = paste0("Kanemaru CMsubtype scores based on t-value marker genes_all CM sub.jpeg"), width=1500 , height =1500,quality = 100,res = 300)
  print(hm)
  dev.off()
  svg(filename = paste0("Kanemaru CMsubtype scores based on t-value marker genes_all CM sub.svg"), width=4,height = 4)
  print(hm)
  dev.off()

#score expression of Reichart CM marker genes in marker genes
  reichart_CM_subtype_markers<-as.data.frame(read_xlsx("S6_CM_Marker_Genes_Cell_States.xlsx"))
  CMsub_states<- unique(sub("_.*", "", colnames(reichart_CM_subtype_markers)))[-1]

  #use logFC as weight for mlm
  for (type in CMsub_states) {
    vec<-paste0(type,c("_gene_symbol","_log2_fold-change"))
    df<-reichart_CM_subtype_markers[,vec]
    colnames(df)<-c("gene_symbol","log2_fold-change")
    df$CM_cs<-paste0("Reichart_",type)
    if (type==CMsub_states[1]) {
      df_combined<-df
    }else{
      df_combined<-rbind(df_combined,df)
    }
  }
  df_combined$FC<-2^df_combined$`log2_fold-change`

  # Run mlm
  sample_acts <- decoupleR::run_mlm(mat = all.markers_full,
                                    net = df_combined,
                                    .source = 'CM_cs',
                                    .target = 'gene_symbol',
                                    .mor = 'FC',
                                    minsize = 1)
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)

  # Color scale
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)

  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05,2, length.out = floor(100 / 2)))
  #add pval
  sample_acts <- sample_acts %>%
    mutate(sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  # Transform to wide format for annotations, similar to `sample_acts_mat`
  annotations <- sample_acts %>%
    tidyr::pivot_wider(id_cols = 'condition',
                       names_from = 'source',
                       values_from = 'sig') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  #plot
  hm<-pheatmap::pheatmap(mat = sample_acts_mat,
                         color = colors.use,
                         cluster_rows = FALSE,
                         cluster_cols = F,
                         border_color = "white",
                         breaks = my_breaks,
                         cellwidth = 20,
                         cellheight = 20,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         display_numbers = annotations,
                         number_color = "black") # Adjust `number_color` as needed
  jpeg(filename = paste0("Reichart CMsubtype scores based on t-value marker genes_all CM sub.jpeg"), width=1500 , height =1500,quality = 100,res = 300)
  print(hm)
  dev.off()
  svg(filename = paste0("Reichart CMsubtype scores based on t-value marker genes_all CM sub.svg"), width=8,height = 4)
  print(hm)
  dev.off()

#======================================================
# DGEA -> plots + GSEA (based on dreamlet result)
#======================================================
#dreamlet result per ct
DEG_list<-qread("per_ct/IVIG_effect.list_list of df per celltype_new.qs")
DEG_list<-DEG_list$CM
DEG_list$gene<-row.names(DEG_list)
#up down factor
DEG_list$up_down<-"up"
DEG_list$up_down[DEG_list$logFC<0]<-"down"
DEG_list$up_down<-factor(DEG_list$up_down,levels = c("up","down"))

#get pseudobulk expression
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg.CM,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_states_label = "subcelltypes.short",
                                     average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)
sc_IVIg_ps<-pseudobulk.QC.filter(sc_IVIg_ps,min.sample = 5)
sc_IVIg_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg_ps@assays$RNA@counts)
sc_IVIg_ps<-NormalizeData_DESeq2(sc_IVIg_ps)
sc_IVIg_ps$three_groups_setup<-as.character(sc_IVIg_ps$condition)
sc_IVIg_ps$three_groups_setup[sc_IVIg_ps$three_groups_setup=="PLACEBO_PRE"|sc_IVIg_ps$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg_ps$three_groups_setup<-factor(sc_IVIg_ps$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))
sc_IVIg_ps<-ScaleData_DESeq2(sc_IVIg_ps,batch1 = "patient_ID",conserve_covariate = "three_groups_setup",regress_covariates = "avg_percent_mt")
sc_IVIg_ps<-FindVariableFeatures_DESeq2(sc_IVIg_ps,nHVG = nrow(sc_IVIg_ps))

#plot only significant
DEG_list_sig<-DEG_list[DEG_list$P.Value<0.05,]

#lable top 15
top10 <- DEG_list_sig[!grepl("^(AC|AP|AL)\\d+", DEG_list_sig$gene),] %>% group_by(up_down) %>% top_n(-15,P.Value) %>% pull(gene)

#plot heatmap
FeatureHeatmap(sc_IVIg_ps,group.by = "condition",features = DEG_list_sig$gene,slot = "scale.data",assay = "DESeq2",group_palcolor = list(col.ivig.condi),
             height = 4.5,width = 1.5,cluster_rows = T,feature_split = DEG_list_sig$up_down,features_label = top10,row_title = c("",""),
             ht_params = list(show_row_dend = FALSE))
ggsave(filename = "DGEA_CM_dreamlet small heatmap.jpeg", width=10 , height = 8)
ggsave(filename = "DGEA_CM_dreamlet small heatmap.svg", width=10 , height = 8)

#enrichment analysis with mlm
# Extract t-values per gene
deg <- DEG_list %>%
  dplyr::select(t) %>%
  dplyr::filter(!is.na(t)) %>%
  as.matrix()

#KEGG as ref
  pw.list.KEGG <- msigdbr("Homo sapiens", "C2","KEGG")
  pw.list.KEGG <- split(pw.list.KEGG$human_gene_symbol, pw.list.KEGG$gs_name)
  pw.list.KEGG <- lapply(pw.list.KEGG, unique)
  KEGG_disease_list <- unlist(read_csv("references/KEGG disease list.txt",col_names = FALSE))
  names(KEGG_disease_list)<-NULL
  KEGG_disease_list <- KEGG_disease_list %>%
    gsub('[0-9]+', '', .) %>%            # Remove numbers
    gsub("\\s*\\([^\\)]*\\)\\s*", " ", .) %>% # Remove text in parentheses
    trimws() %>%                          # Trim leading and trailing whitespace
    toupper() %>%                         # Convert to uppercase
    gsub("[- ]", "_", .) %>%              # Replace hyphens and spaces with underscores
    gsub("_{2,}", "_", .) %>%             # Replace multiple underscores with a single one
    paste0("KEGG_", .)                    # Add "KEGG_" prefix
  pw.list.KEGG<-pw.list.KEGG[!names(pw.list.KEGG) %in% KEGG_disease_list]
  pw.list.KEGG<-pw.list.KEGG[!str_detect(names(pw.list.KEGG),"INFECTION")]
  pw.list.KEGG<-pw.list.KEGG[!str_detect(names(pw.list.KEGG),"DISEASE")]
  pw.list.KEGG<-pw.list.KEGG[!str_detect(names(pw.list.KEGG),"CARDIOMYOPATHY")]
  pw.list.KEGG<-pw.list.KEGG[!str_detect(names(pw.list.KEGG),"DEPRESSION")]
  net<-melt(pw.list.KEGG)
  net<-net[!duplicated(net),]
  net$mor<-1

  # Run mlm
  contrast_acts <- decoupleR::run_mlm(mat = deg,
                                    net = net,
                                    .source = 'L1',
                                    .target = 'value',
                                    minsize = 5)
  contrast_acts<-arrange(contrast_acts,p_value)
  contrast_acts<-contrast_acts[contrast_acts$p_value<0.1,]
  contrast_acts<-arrange(contrast_acts,score)

  #dotplot
  ggplot(contrast_acts, aes(x=stats::reorder(source, score), y="KEGG mlm",size=-log10(p_value),fill=score)) +
    geom_point(shape = 21, stroke = 0.5,) +
    scale_size(range = c(2,6)) +
    scale_x_discrete(position = 'top')+
    coord_flip()+
    theme_classic()+
    scale_fill_gradientn(colours = special.ramp(contrast_acts$score))
  ggsave(filename = "mlm of t-value from dreamlet on CM_KEGG _dp.jpeg", width=6, height = 4)
  ggsave(filename = "mlm of t-value from dreamlet on CM_KEGG _dp.svg", width=6 , height = 4)

  #expression score dotplot in the snRNA data
  sig_terms<-pw.list.KEGG[names(pw.list.KEGG) %in% contrast_acts$source]
  sc_IVIg.CM<-AddModuleScore(sc_IVIg.CM,features = sig_terms,name = names(sig_terms))
  colnames(sc_IVIg.CM@meta.data)[(length(colnames(sc_IVIg.CM@meta.data))-(length(names(sig_terms))-1)):length(colnames(sc_IVIg.CM@meta.data))]<-names(sig_terms)
  DotPlot(sc_IVIg.CM,group.by = "subcelltypes.short",features = contrast_acts$source)+coord_flip()+xlab("")+ylab("")+ theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(filename = "AddModuleScore expression of sig_KEGG_terms _dp.jpeg", width=7.2, height = 5)
  ggsave(filename = "AddModuleScore expression of sig_KEGG_terms _dp.svg", width=7.2 , height = 5)
