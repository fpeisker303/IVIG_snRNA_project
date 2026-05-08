#IVIG project
#pathway centric analysis in fibroblast

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")
library(rstatix)

#=========================================================================================
#load NABA gene set
#=========================================================================================
  # Process human NABA matrisome genes from #http://matrisomeproject.mit.edu/
  naba_ECM <- as.data.frame(read_excel("Hs_Matrisome_Masterlist_Naba et al_2012.xlsx"))
  naba_ECM$`Matrisome Division`=gsub(pattern = " ",replacement = "_",x = naba_ECM$`Matrisome Division`)
  naba_ECM$`Matrisome Division`=gsub(pattern = "-",replacement = "_",x = naba_ECM$`Matrisome Division`)
  naba_ECM$`Matrisome Category`=gsub(pattern = " ",replacement = "_",x = naba_ECM$`Matrisome Category`)
  naba_ECM$`Matrisome Category`=gsub(pattern = "-",replacement = "_",x = naba_ECM$`Matrisome Category`)
  naba_ECM.list=NULL
  for (geneset in unique(naba_ECM$`Matrisome Category`)) {
    naba_ECM.list[[geneset]] = naba_ECM$`Gene Symbol`[naba_ECM$`Matrisome Category`==geneset]
  }
  for (geneset in unique(naba_ECM$`Matrisome Division`)) {
    naba_ECM.list[[geneset]] = naba_ECM$`Gene Symbol`[naba_ECM$`Matrisome Division`==geneset]
  }
  #Core_matrisome and Matrisome_associated summarize the other terms and need to excluded here
  net<-melt(naba_ECM.list[1:6])
  net<-net[!duplicated(net),]
  net$mor<-1
  
#=========================================================================================
## in FIB
# score immune terms with mlm in corrected/scaled pseudobulk matrix per sample
# calculate log(ratio) per patient_ID
# compare PLACEO vs IVIG
#=========================================================================================
#load FIB
sc_IVIg.FIB<-qread("subclustering of FB_annotated.qs")

#get pseudobulk expression matrix for scoring
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg.FIB,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",cell_number_cutoff = 5,
                                     average_meta_feature = c("percent.mt","percent.rb","doublet_score"),add_ct_cs_faction = T)
#filter out genes that are only found very low expressed in a very few samples
sc_IVIg_ps<-pseudobulk.QC.filter(sc_IVIg_ps,min.percent.sample =  0.5) 
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

df<-dcast(sample_acts,source ~ condition,value.var = "score")
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
  df_tmp<-apply(df_tmp, 1, function(x){return(x[2]-x[1])})
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
#stat_test <- stat_test[stat_test$p<0.05,] #keep all terms, indicate significance on the plot
stat_test$term<-factor(stat_test$term,levels = stat_test$term)

stat_test <-  add_y_position(stat_test,step.increase = 0.15,fun = "max")
stat_test$x<-1:nrow(stat_test)
stat_test$xmin<-stat_test$x-0.2
stat_test$xmax<-stat_test$x+0.2

#optimize y position
max_y <- df_differences %>% 
  group_by(term) %>%
  summarise(max_value = max(value, na.rm = TRUE)) %>% as.data.frame()
max_y$y.position<-max_y$max_value + max_y$max_value*0.03
max_y$y.position[max_y$max_value < 5]<-max_y$max_value[max_y$max_value < 5]+0.2
rownames(max_y)<-max_y$term
#override
stat_test$y.position <- max_y[as.character(stat_test$term),]$y.position

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
  ylab("Expression score difference per group")+
  xlab("")+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_classic()+
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=10))+
  stat_pvalue_manual(stat_test,label = "p.adj", tip.length = 0.01,coord.flip = T,size = 3)+coord_flip()

ggsave(filename = paste0("mlm score diff ps based per sample for NABA ECM.jpeg"), width=6 , height = 5,limitsize = F)
ggsave(filename = paste0("mlm score diff ps based per sample for NABA ECM.svg"), width=6 , height = 5,limitsize = F)