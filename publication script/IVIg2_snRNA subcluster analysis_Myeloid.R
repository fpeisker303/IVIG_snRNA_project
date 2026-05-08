## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#subcluster analysis for myeloid cells

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#load
sc_IVIg.Myel<-qread("./subclustering of Myel_annotated_f5.qs")

#add the extended major ct annotation of myeloid
sc_IVIg.Myel$celltypes.short.ext<-as.character(sc_IVIg.Myel$subcelltypes.short)
sc_IVIg.Myel$celltypes.short.ext[str_detect(sc_IVIg.Myel$celltypes.short.ext,"MP_")]<-"Macrophages"
sc_IVIg.Myel$celltypes.short.ext[sc_IVIg.Myel$celltypes.short.ext %in% c("MO_CD16", "MO_VCAN" )]<-"Monocytes"
sc_IVIg.Myel$celltypes.short.ext<-factor(sc_IVIg.Myel$celltypes.short.ext,levels = c("cDC","Monocytes","Macrophages","MAST"))

#=========================================================================================
# gene expression analysis
#=========================================================================================
#FC receptors
#FcyRIA, FcyRIIA, FcyRIIC, FcyRIIIA (CD16), FcyRIIIB (CD16b), FcyRIIB, FcRn, DC-SIGN
#FCGR1A, FCGR2A, FCGR2C, FCGR3A, FCGR3B, FCGR2B, FCGRT, CD209
FcyR_genes<-c("FCGR1A", "FCGR2A", "FCGR3A", 'FCGR3B', "FCGR2B", "FCGRT", "CD209") #"FCGR2C" not in data
sc_IVIg.Myel<-ScaleData(sc_IVIg.Myel,features = FcyR_genes)

DotPlot(sc_IVIg.Myel,features = FcyR_genes,group.by = "celltypes.short.ext")+coord_flip()  +theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0("Fc receptor expressions in major ct_dp.jpeg") , width=4.5, height = 3.2)
ggsave(filename = paste0("Fc receptor expressions in major ct_dp.svg") , width=4.5, height = 3.2)

#======================================================
# DGEA -> plots + GSEA (based on dreamlet result) 
#======================================================
#dreamlet result per ct
DEG_list<-qread("IVIG_effect.list_list of df per celltype_new.qs")
DEG_list<-DEG_list$Myeloid
DEG_list$gene<-row.names(DEG_list)
#up down factor
DEG_list$up_down<-"up"
DEG_list$up_down[DEG_list$logFC<0]<-"down"
DEG_list$up_down<-factor(DEG_list$up_down,levels = c("up","down"))

#get pseudobulk expression
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg.Myel,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_states_label = "subcelltypes.short",
                                     average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)
sc_IVIg_ps<-pseudobulk.QC.filter(sc_IVIg_ps,min.sample = 4) 
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
top10 <- DEG_list_sig[!grepl("^(AC|AP|AL)\\d+", DEG_list_sig$gene),] %>% group_by(up_down) %>% top_n(-10,P.Value) %>% pull(gene)

#plot heatmap
FeatureHeatmap(sc_IVIg_ps,group.by = "condition",features = DEG_list_sig$gene,slot = "scale.data",assay = "DESeq2",group_palcolor = list(col.ivig.condi),
               height = 4.5,width = 1.5,cluster_rows = T,feature_split = DEG_list_sig$up_down,features_label = top10,row_title = c("",""),
               ht_params = list(show_row_dend = FALSE))
ggsave(filename = "DGEA_Myel_dreamlet small heatmap.jpeg", width=10 , height = 8)
ggsave(filename = "DGEA_Myel_dreamlet small heatmap.svg", width=10 , height = 8)

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
ggsave(filename = "mlm of t-value from dreamlet on Myel_KEGG _dp.jpeg", width=6.1 , height = 4)
ggsave(filename = "mlm of t-value from dreamlet on Myel_KEGG _dp.svg", width=6.1 , height = 4)

#expression score dotplot in the snRNA data
sig_terms<-pw.list.KEGG[names(pw.list.KEGG) %in% contrast_acts$source]
sc_IVIg.Myel<-AddModuleScore(sc_IVIg.Myel,features = sig_terms,name = names(sig_terms))
colnames(sc_IVIg.Myel@meta.data)[(length(colnames(sc_IVIg.Myel@meta.data))-19):length(colnames(sc_IVIg.Myel@meta.data))]<-names(sig_terms)
DotPlot(sc_IVIg.Myel,group.by = "celltypes.short.ext",features = contrast_acts$source)+coord_flip()+xlab("")+ylab("")+ theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = "AddModuleScore expression of sig_KEGG_terms _dp.jpeg", width=6.9 , height = 4.3)
ggsave(filename = "AddModuleScore expression of sig_KEGG_terms _dp.svg", width=6.9 , height = 4.3)
