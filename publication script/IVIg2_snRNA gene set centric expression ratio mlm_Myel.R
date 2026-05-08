#IVIG project
#pathway centric analysis in myeloid

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")
library(rstatix)

#=========================================================================================
#load gene set db and select immunologic functional processes
#=========================================================================================
pw.list.KEGG <- msigdbr("Homo sapiens", "C2","KEGG")
pw.list.GOBP <- msigdbr("Homo sapiens", "C5","GO:BP")

#read immune terms from KEGG 
KEGG_immune<-as.data.frame(read.delim("immune related KEGG 09151.txt",header = F))
KEGG_immune<-toupper(KEGG_immune$V1)
KEGG_immune<-paste0("KEGG_",gsub(" ","_",KEGG_immune))
KEGG_immune<-gsub("-","_",KEGG_immune)
#remove terms with Lymphoid cell functions
KEGG_immune<-KEGG_immune[!str_detect(pattern = "T_CELL|B_CELL|NATURAL_KILLER|TH17|TH1|IGA_PRODUCTION",string = KEGG_immune)]
# Convert to list of pathways
pw.list.KEGG <- split(pw.list.KEGG$human_gene_symbol, pw.list.KEGG$gs_name)
pw.list.KEGG <- lapply(pw.list.KEGG, unique)
#select immune terms from KEGG 
pw.list.KEGG<-pw.list.KEGG[KEGG_immune][!is.na(names(pw.list.KEGG[KEGG_immune]))]

# look for immune terms from GOPB
terms_GOBP<-unique(pw.list.GOBP$gs_name[str_detect(pattern ="IMMUNE",string =  pw.list.GOBP$gs_name)])
terms_GOBP<-c(terms_GOBP,unique(pw.list.GOBP$gs_name[str_detect(pattern ="CYTOKINE",string =  pw.list.GOBP$gs_name)]))
terms_GOBP<-unique(c(terms_GOBP,unique(pw.list.GOBP$gs_name[str_detect(pattern ="CHEMOKINE",string =  pw.list.GOBP$gs_name)])))
#remove terms with Lymphoid cell functions
terms_GOBP<-terms_GOBP[!str_detect(pattern = "T_CELL|T_HELPER|NATURAL_KILLER|CD4|B_CELL|SOMATIC_RECOMBINATION",string = terms_GOBP)]
#prepare geneset lists and subset
pw.list.GOBP <- split(pw.list.GOBP$human_gene_symbol, pw.list.GOBP$gs_name)
pw.list.GOBP <- lapply(pw.list.GOBP, unique)
#select immune terms from GOPB 
pw.list.GOBP<-pw.list.GOBP[terms_GOBP]

#combine both lists
pw.list_immune<-c(pw.list.GOBP,pw.list.KEGG)

#remove terms with cytokinesis, as this has nothing specific to do with cytokines
pw.list_immune<-pw.list_immune[!str_detect(names(pw.list_immune),"CYTOKINESIS")]

#input for mlm
net<-melt(pw.list_immune)
net<-net[!duplicated(net),]
net$mor<-1

#some terms have high colinearity
colinear_res<-decoupleR::check_corr(network =net,.source = 'L1',.target = 'value')
#remove one term for correlation >0.7
colinear_res <- colinear_res[colinear_res$correlation>0.7,]
#rm terms with high correlation
net<-net[!net$L1%in%colinear_res$L1.2,]

#=========================================================================================
## in Myeloid
# score immune terms with mlm in corrected/scaled pseudobulk matrix per sample
# calculate log(ratio) per patient_ID
# compare PLACEBO vs IVIG
#=========================================================================================
#load Myeloid
sc_IVIg.Myel<-qread("subclustering of Myel_annotated_f5.qs")

#get pseudobulk expression matrix for scoring
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg.Myel,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_number_cutoff = 10,
                                     average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)
#filter out genes that are only found very low expressed in a very few samples
sc_IVIg_ps<-pseudobulk.QC.filter(sc_IVIg_ps,min.sample = 3) 
#Create a new assay to process the data with DESeq2 functions
sc_IVIg_ps[["DESeq2"]]<-CreateAssayObject(counts = sc_IVIg_ps@assays$RNA@counts)
sc_IVIg_ps<-NormalizeData_DESeq2(sc_IVIg_ps)
#correct for avg_percent_mt and patient_ID with the 3 groups design
sc_IVIg_ps$three_groups_setup<-as.character(sc_IVIg_ps$condition)
sc_IVIg_ps$three_groups_setup[sc_IVIg_ps$three_groups_setup=="PLACEBO_PRE"|sc_IVIg_ps$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg_ps$three_groups_setup<-factor(sc_IVIg_ps$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))
sc_IVIg_ps<-ScaleData_DESeq2(sc_IVIg_ps,batch1 = "patient_ID",conserve_covariate = "three_groups_setup",regress_covariates = "avg_percent_mt")
meta.data<-sc_IVIg_ps@meta.data

# Run mlm
sample_acts <- decoupleR::run_mlm(mat = sc_IVIg_ps@assays$DESeq2$scale.data, 
                                  net = net, 
                                  .source = 'L1', 
                                  .target = 'value',
                                  minsize = 5)

#change enrichment score to positive score
sample_acts$score_pos<-sqrt(10^sample_acts$score)

df<-dcast(sample_acts,source ~ condition,value.var = "score_pos")
rownames(df)<-paste0(df$source)
df$source<-NULL

#calculate ratio of the gsva score per paired sample
#get samples with pairs
samples<-unique(rownames(meta.data))
samples<-gsub("_PRE","",samples)
samples<-gsub("_POST","",samples)
meta.data$patient_ID_2<-samples
samples_with_pairs<-names(table(samples))[  table(samples)>1]
meta.data<-meta.data[meta.data$patient_ID_2 %in% samples_with_pairs,]

#reduce to pairs
df<-df[,as.character(rownames(meta.data))]

#calculate ratio per patient
df_differences<-NULL
for (patient in unique(meta.data$patient_ID_2)) {
  df_tmp<-df[,rownames(meta.data[meta.data$patient_ID_2==patient,])]
  df_tmp<-apply(df_tmp, 1, function(x){return(log(x[2]/x[1]))})
  new.df<-data.frame("patient_ID"=patient,"term" = names(df_tmp),"value" =df_tmp)
  if (patient==unique(meta.data$patient_ID_2)[1]) { df_differences<-new.df }else{df_differences<-rbind(df_differences,new.df)}
}
rownames(df_differences)<-NULL
#get treatment groups from meta.data
df_differences$study_group <- meta.data$study_group[match(df_differences$patient_ID, meta.data$patient_ID)]

# Statistical test
# Perform pairwise comparisons using Wilcoxon tests
stat_test <- df_differences  %>%
  group_by(term) %>%
  pairwise_wilcox_test(value ~ study_group, p.adjust.method = "BH")  %>%
  add_significance("p.adj")
stat_test<-arrange(stat_test,p)
stat_test <- stat_test[stat_test$p.adj<0.05,] #filter for significance
stat_test$term<-factor(stat_test$term,levels = stat_test$term)

stat_test <-  add_y_position(stat_test,step.increase = 0.15,fun = "max")
stat_test$x<-1:nrow(stat_test)
stat_test$xmin<-stat_test$x-0.2
stat_test$xmax<-stat_test$x+0.2

#add column with only significant pvalues
stat_test$p.adj.signif.num<-stat_test$p.adj
stat_test$p.adj.signif.num[stat_test$p.adj.signif.num > 0.05]<-"ns"
stat_test$p.adj.signif.NA<-stat_test$p.adj
stat_test$p.adj.signif.NA[stat_test$p.adj.signif.NA > 0.05]<-NA

#get significant terms
ratios<-df_differences[df_differences$term %in%stat_test$term,]
ratios$term<-factor(ratios$term,levels=stat_test$term)
ratios<-arrange(ratios,term)

#plot significant terms    
ggplot(ratios, aes(x=term, y=value)) + 
  geom_boxplot(aes(fill=study_group),position=position_dodge2(padding=0.1), pch=21, show.legend = T,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.8),pch=21, aes(fill=study_group), show.legend = T)+
  scale_fill_manual(values =col.ivig.condi[c(2,4)])+
  ylab("Expression score logFC")+
  xlab("")+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=10))+
  stat_pvalue_manual(stat_test,label = "p.adj.signif", tip.length = 0.01,coord.flip = T,size = 6)+coord_flip()

ggsave(filename = paste0("mlm score ps based per sample for immune terms_myeloid.jpeg"), width=9 , height = 3,limitsize = F)
ggsave(filename = paste0("mlm score ps based per sample for immune terms_myeloid.svg"), width=9 , height = 3,limitsize = F)