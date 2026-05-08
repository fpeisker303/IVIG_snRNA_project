#requires R4.3 or higher for dreamlet
#IVIg DCM snRNA project
#perform DGEA based on pseudobulk the data

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

library(SingleCellExperiment)
library(dreamlet)
library(muscat)
library(scater)

#load
sc_IVIg<-qread("IVIG_snRNA_harmony_f5.qs")

#=========================================================================================
# DGEA on major celltype. except for LYMEC and ADI, because only very low numbers
#=========================================================================================
#filter 
Idents(sc_IVIg)="celltypes.short"
sc_IVIg<-subset(sc_IVIg,idents=c("LYMEC","ADI","NEURO"),invert=T)
sc_IVIg$celltypes.short<-factor(sc_IVIg$celltypes.short)

sc_IVIg$three_groups_setup<-as.character(sc_IVIg$condition)
sc_IVIg$three_groups_setup[sc_IVIg$three_groups_setup=="PLACEBO_PRE"|sc_IVIg$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg$three_groups_setup<-factor(sc_IVIg$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))


#=========================================================================================
# using dreamlet to perform DGEA on in all cell types 
#=========================================================================================
#https://diseaseneurogenomics.github.io/dreamlet/index.html
sce<-as.SingleCellExperiment(sc_IVIg)

# Create pseudobulk data by specifying cluster_id and sample_id
# Count data for each cell type is then stored in the `assay` field
# assay: entry in assayNames(sce) storing raw counts
# cluster_id: variable in colData(sce) indicating cell clusters
# sample_id: variable in colData(sce) indicating sample id for aggregating cells
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "celltypes.short",
                            sample_id = "orig.ident",
                            verbose = FALSE
)

# one 'assay' per cell type
assayNames(pb)

#add meta data to mt
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_states_label = "subcelltypes.short",
                                         average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)

pb$avg_percent_mt<-sc_IVIg_ps$avg_percent_mt

# Normalize and apply voom/voomWithDreamWeights
formula <- ~0 + three_groups_setup + avg_percent_mt + (1|patient_ID)
param = SnowParam(14, "SOCK", progressbar=TRUE,exportglobals = FALSE)
res.proc <- processAssays(pb, formula , min.count = 5,BPPARAM = param)

qsave(res.proc,"dreamlet_res.proc_3groups.qs")
#res.proc<-qread("./dreamlet_res.proc_3groups.qs")

# Differential expression analysis within each assay,
  # evaluated on the voom normalized data
  res.dl <- dreamlet(res.proc,formula = formula,BPPARAM = param,contrasts = c(IVIG_effect = "three_groups_setupIVIG_POST - three_groups_setupPLACEBO_POST"))
  
#get results per cell type
  IVIG_effect.list<-NULL
  #plot through celltypes for heatmap
  for (celltype in unique(sc_IVIg$celltypes.short)) {
    message(celltype)
    #adjust folder
    setwd(paste0("./per_ct"))
    dir.create(celltype)
    setwd(paste0("./",celltype))
    
    #subset
    Idents(sc_IVIg)<-"celltypes.short"
    sc_IVIg_sub<-subset(sc_IVIg,ident=celltype)
    
    #create the pseudobulk + filterd and get meta.data ref
    sc_IVIg_sub_ps<-create_pseudobulk_object(sc_IVIg_sub,sample_ID = "orig.ident",average_meta_feature = c("percent.mt","percent.rb","doublet_score"))
    sc_IVIg_sub_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg_sub_ps@assays$RNA@counts)
    sc_IVIg_sub_ps<-NormalizeData_DESeq2(sc_IVIg_sub_ps)  
    
    #filter out genes that are only found very low expressed in a very few samples
    sc_IVIg_sub_ps<-pseudobulk.QC.filter(sc_IVIg_sub_ps,min.sample = 5) 
    sc_IVIg_sub_ps<-NormalizeData_DESeq2(sc_IVIg_sub_ps)  
    
    ## get the result for IVIG_effect
    res_IVIG_effect <- topTable(res.dl[[celltype]],number = dim(res.dl[[celltype]])[1], coef = "IVIG_effect")
    res_IVIG_effect_full<-res_IVIG_effect
    res_IVIG_effect<-res_IVIG_effect[res_IVIG_effect$P.Value<0.05,]
    res_IVIG_effect<-arrange(res_IVIG_effect,P.Value)
    res_IVIG_effect<-arrange(res_IVIG_effect,-logFC)
    res_IVIG_effect$gene<-rownames(res_IVIG_effect)
    res_IVIG_effect$up_down<-res_IVIG_effect$logFC<0
    
    #get the top 20 down and up to label in the heatmap (exclude lnRNA from the labels)
    res_IVIG_effect_up<-res_IVIG_effect[res_IVIG_effect$logFC>0 & res_IVIG_effect$P.Value<0.05,]
    res_IVIG_effect_up<-arrange(res_IVIG_effect_up,P.Value)
    res_IVIG_effect_up <- res_IVIG_effect_up[!grepl("^(AC|AP|AL|LINC)\\d+", res_IVIG_effect_up$gene), ] #filter lnRNA
    res_IVIG_effect_up<-res_IVIG_effect_up[1:20,]
    res_IVIG_effect_down<-res_IVIG_effect[res_IVIG_effect$logFC<0 & res_IVIG_effect$P.Value<0.05,]
    res_IVIG_effect_down<-arrange(res_IVIG_effect_down,P.Value)
    res_IVIG_effect_down <- res_IVIG_effect_down[!grepl("^(AC|AP|AL|LINC)\\d+", res_IVIG_effect_down$gene), ] #filter lnRNA
    res_IVIG_effect_down<-res_IVIG_effect_down[1:20,]
    res_IVIG_effect_top<-rbind(res_IVIG_effect_up,res_IVIG_effect_down)
    
    #save result for GSEA
    IVIG_effect.list[[celltype]]<-res_IVIG_effect_full
  
  }
qsave(IVIG_effect.list,"IVIG_effect.list_list of df per celltype_new.qs")  
