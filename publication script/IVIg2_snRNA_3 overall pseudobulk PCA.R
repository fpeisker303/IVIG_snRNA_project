#IVIG project
#overall pseudobulk PCA

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#load data
sc_IVIg<-qread("IVIG_snRNA_harmony_f5.qs")

#=========================================================================================
# correction for patient_ID (3 group design) and percent_mt
# in addition correcting for the endoEC percent
#=========================================================================================
#there is a strong influence due to endoEC proportion, especially higher in IVIG_POST
#since this is clearly a sampling bias, trying regression for the pct of EndoEC

#load ps and add 3 group design
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_states_label = "subcelltypes.short",
                                     average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)
sc_IVIg_ps<-pseudobulk.QC.filter(sc_IVIg_ps,min.sample = 5) 
sc_IVIg_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg_ps@assays$RNA@counts)
sc_IVIg_ps<-NormalizeData_DESeq2(sc_IVIg_ps)

#add 3 groups design
sc_IVIg_ps$three_groups_setup<-as.character(sc_IVIg_ps$condition)
sc_IVIg_ps$three_groups_setup[sc_IVIg_ps$three_groups_setup=="PLACEBO_PRE"|sc_IVIg_ps$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg_ps$three_groups_setup<-factor(sc_IVIg_ps$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))

#correction
sc_IVIg_ps_sub<-ScaleData_DESeq2(sc_IVIg_ps,batch1 = "patient_ID",conserve_covariate = "three_groups_setup",regress_covariates = c("avg_percent_mt","ct_EndoEC_pct"))

#PCA
sc_IVIg_ps_sub<-FindVariableFeatures_DESeq2(sc_IVIg_ps_sub,nHVG = nrow(sc_IVIg_ps_sub))
sc_IVIg_ps_sub<-RunPCA_DESeq2(sc_IVIg_ps_sub)

#Plot
CellDimPlot(sc_IVIg_ps_sub,reduction = "DESeq2_pca",group.by = "condition",pt.size = 3,palcolor = col.ivig.condi,add_mark = T,mark_type = "ellipse")
ggsave(filename = paste0("DESeq2_pca_condition_ellipse with correction.jpeg"), width=8, height = 7,limitsize = F)
ggsave(filename = paste0("DESeq2_pca_condition_ellipse with correction.svg"), width=8, height = 7,limitsize = F)
