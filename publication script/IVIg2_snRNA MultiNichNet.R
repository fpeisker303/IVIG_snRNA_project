#IVIG project
#cell cell communication analysis with multinichenet v2.0
#important alteration from the multinichenet default pipeline!!! 
#here: differential expression result of dreamlet used instead of muscat!!!

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

library(SingleCellExperiment)
library(nichenetr)
library(multinichenetr)

#=============================================================
# load the nichenet LR model
#=============================================================
organism = "human"
options(timeout = 120)
if(organism == "human"){
  
  lr_network_all = 
    readRDS("references/lr_network_human_allInfo_30112033.rds") %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS("references/ligand_target_matrix_nsga2r_final.rds")
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

#=============================================================
# load and prepare data. and input parameter
#=============================================================
#load snRNA data
sc_IVIg<-qread("IVIG_snRNA_harmony_f5.qs")

#add 3 groups design
sc_IVIg$three_groups_setup<-as.character(sc_IVIg$condition)
sc_IVIg$three_groups_setup[sc_IVIg$three_groups_setup=="PLACEBO_PRE"|sc_IVIg$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg$three_groups_setup<-factor(sc_IVIg$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))

#create a pseudobulk to get the avg_percent_mt in order to add this to the snRNA data as covariate
sc_IVIg_ps<-create_pseudobulk_object(sc_IVIg,sample_ID = "orig.ident",average_meta_feature = c("percent.mt","percent.rb"))
df<-left_join(FetchData(sc_IVIg,vars = c("orig.ident")),FetchData(sc_IVIg_ps,vars = c("orig_ident","avg_percent_mt")),by = c("orig.ident"="orig_ident"))
rownames(df)<-rownames(FetchData(sc_IVIg,vars = c("orig.ident")))
sc_IVIg<-AddMetaData(sc_IVIg,df[2],col.name = "avg_percent_mt")

#convert seurat to sce and check gene names
sce = Seurat::as.SingleCellExperiment(sc_IVIg, assay = "RNA")
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
rm(sc_IVIg,sc_IVIg_ps)
gc()

#add info for meta data columns to use
sample_id = "orig.ident"
group_id = "three_groups_setup"
celltype_id = "celltypes.short"

#covariates will just be included in the DE GLM model
#batches will be included in the DE GLM model AND normalized pseudobulk expression values will be corrected for the batch effects
covariates = c("patient_ID","avg_percent_mt")
batches = NA
#batches is handled in the DE also simply as a covariate in the design, therefore patient_ID is used here as covariate
#first run MNN DE as implemented, potentially replace the result with the one from dreamlet!

#define the contrast for DGEA
contrasts_oi = c("'IVIG_POST-PLACEBO_POST'")
#For downstream visualizations and linking contrasts to their main condition, we also need to run the following
contrast_tbl = tibble(contrast = c("IVIG_POST-PLACEBO_POST","PLACEBO_POST-IVIG_POST"),
                      group = c("IVIG_POST", "PLACEBO_POST"))

#Define the sender and receiver cell types of interest
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

#=============================================================
#Cell-type filtering: determine which cell types are sufficiently present
#=============================================================
#min cells for pseudobulk
min_cells = 10
#check abundance
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
)
#plot
abundance_info$abund_plot_sample

#Cell type filtering based on cell type abundance information
sample_group_celltype_df = abundance_info$abundance_data %>% 
  filter(n > min_cells) %>% 
  ungroup() %>% 
  distinct(sample_id, group_id) %>% 
  cross_join(
    abundance_info$abundance_data %>% 
      ungroup() %>% 
      distinct(celltype_id)
  ) %>% 
  arrange(sample_id)

abundance_df = sample_group_celltype_df %>% left_join(
  abundance_info$abundance_data %>% ungroup()
)

abundance_df$n[is.na(abundance_df$n)] = 0
abundance_df$keep[is.na(abundance_df$keep)] = FALSE
abundance_df_summarized = abundance_df %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

# absent_celltypes = abundance_df_summarized %>% #DOES NOT WORK
#   filter(samples_present < 2) %>% 
#   group_by(celltype_id) %>% 
#   count() %>% 
#   filter(n == total_nr_conditions) %>% 
#   pull(celltype_id)

print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## character(0)

print("absent celltypes:")
## [1] "absent celltypes:"
# print(absent_celltypes)
## character(0)

# ## filter #DOES NOT WORK, also not needed here
# analyse_condition_specific_celltypes = FALSE
# if(analyse_condition_specific_celltypes == TRUE){
#   senders_oi = senders_oi %>% setdiff(absent_celltypes)
#   receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
# } else {
#   senders_oi = senders_oi %>% 
#     setdiff(union(absent_celltypes, condition_specific_celltypes))
#   receivers_oi = receivers_oi %>% 
#     setdiff(union(absent_celltypes, condition_specific_celltypes))
# }

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]
#=============================================================
#Gene filtering: determine which genes are sufficiently expressed in each present cell type
#=============================================================
# genes should be expressed in at least 2 samples if the group with lowest nr. of samples has 4 samples like this dataset.
min_sample_prop = 0.50

#genes should show non-zero expression values in at least 5% of cells in a sample.
fraction_cutoff = 0.05

#calculate the information required for gene filtering 
frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

#only keep genes that are expressed by at least one cell type
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]


#=============================================================
#Pseudobulk expression calculation: determine and normalize per-sample pseudobulk expression levels for each expressed gene in each present cell type
#=============================================================
abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)
empirical_pval = FALSE

#=============================================================
#Differential expression (DE) analysis: determine which genes are differentially expressed
#=============================================================
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)

#plot pval dist
DE_info$hist_pvals

#In case these p-value distributions look irregular, you can estimate empirical p-values as we will demonstrate in another vignette.
empirical_pval = FALSE
#get celltype de info
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

#check celltype_de output as the DGEA result used for the downstream analysis!
celltype_de

#OPTIONAL: remove the DEG of EndoEC since they are likley overrepresented in the data and in some sample especially
celltype_de<-celltype_de[!celltype_de$cluster_id=="EndoEC",]

#compare to dreamlet results
DEG_list<-qread("/per_ct/IVIG_effect.list_list of df per celltype_new.qs")
#merge into one dataframe
#replace the muscat result by the results from dreamlet
for (celltype in names(DEG_list)) {
  DEG_df_new<-DEG_list[[celltype]]
  DEG_df_new$cluster_id<-celltype
  DEG_df_new$gene<-rownames(DEG_df_new)
  #later on only gene, cluster_id, logFC, p_val, p_adj, contrast are required
  colnames(DEG_df_new)<-c("logFC" ,"AveExpr","t" , "p_val","p_adj","B" , "z.std" ,"cluster_id" ,"gene")
  rownames(DEG_df_new)<-NULL
  #correct cell names in dreamlet result
  DEG_df_new$gene<-convert_alias_to_symbols(DEG_df_new$gene,organism = "human")
  #reduce dreamlet result to genes not filtered by the MNN pipeline beforehand
  MNN_genes<-celltype_de$gene[celltype_de$cluster_id==celltype]
  DEG_df_new<-DEG_df_new[DEG_df_new$gene %in% MNN_genes,]
  if (celltype==names(DEG_list)[1]) { 
    DEG_df<-DEG_df_new
    }else{
    DEG_df<-rbind(DEG_df,DEG_df_new)
  }
}
#this contrast shows upregulation in IVIG
DEG_df$contrast<-"IVIG_POST-PLACEBO_POST"
#add the inverse contrast by changeing the sign to also get what is PLACEBO specific in the later downstream analysis
DEG_df2<-DEG_df
DEG_df2$logFC<-DEG_df2$logFC*-1
DEG_df2$contrast<-"PLACEBO_POST-IVIG_POST"
DEG_df<-rbind(DEG_df,DEG_df2)
celltype_de<-as_tibble(DEG_df)

#reduce sender and receiver oi to celltypes with sig. DEG
senders_oi<-senders_oi[senders_oi %in% unique(celltype_de$cluster_id)]
receivers_oi<-receivers_oi[receivers_oi %in% unique(celltype_de$cluster_id)]
senders_oi<-factor(senders_oi)
receivers_oi<-factor(receivers_oi)

#Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

#=============================================================
#Ligand activity prediction:
#use the DE analysis output to predict the activity of ligands in receiver cell types and infer their potential target genes
#=============================================================
#Assess geneset_oi-vs-background ratios for different DE output tresholds prior to the NicheNet ligand activity analysis
#To determine the genesets of interest based on DE output, we need to define some logFC and/or p-value thresholds per cell type/contrast combination.
#In general, we recommend inspecting the nr. of DE genes for all cell types based on the default thresholds and adapting accordingly.
logFC_threshold = 0.50
p_val_threshold = 0.05
#By default, we will apply the p-value cutoff on the normal p-values, and not on the p-values corrected for multiple testing. 
#This choice was made because most multi-sample single-cell transcriptomics datasets have just a few samples per group and we might have a lack of statistical power due to pseudobulking. 
p_val_adj = FALSE 

#check now if all cell type / contrast combinations, all geneset/background ratio's are within the recommended range
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

#Perform the ligand activity analysis and ligand-target inference
#select which top n of the predicted target genes will be considered (here: top 250 targets per ligand).
#We recommend users to test other settings in case they would be interested in exploring fewer, but more confident target genes, or vice versa.
top_n_target = 250

#set up parameter for parallel
verbose = TRUE
cores_system = 12
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

#run the ligand activity prediction
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))

#=============================================================
#Prioritization: rank cell-cell communication patterns through multi-criteria prioritization
#=============================================================
#Here we will choose for setting ligand_activity_down = FALSE and focus specifically on upregulating ligands.
ligand_activity_down = F

#We will combine these prioritization criteria in a single aggregated prioritization score. 
#In the default setting, we will weigh each of these criteria equally (scenario = "regular"). 
#The setting scenario = "lower_DE" halves the weight for DE criteria and doubles the weight for ligand activity. 
#The setting scenario = "no_frac_LR_expr" ignores the criterion "Sufficiently high expression levels of ligand and receptor in many samples of the same group". 
#This may be interesting for users that have data with a limited number of samples and don’t want to penalize interactions if they are not sufficiently expressed in some samples.
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

#Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
)

#Save all the output of MultiNicheNet
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)
qsave(multinichenet_output, paste0("multinichenet_output.qs"))

#=============================================================
# load results and parameter data of MNN 
#=============================================================
multinichenet_output<-qread("./multinichenet_output.qs")

#parameter info
sample_id = "orig.ident"
group_id = "three_groups_setup"
celltype_id = "celltypes.short"
covariates = c("patient_ID","avg_percent_mt")
batches = NA
contrasts_oi = c("'IVIG_POST-PLACEBO_POST'")
contrast_tbl = tibble(contrast = c("IVIG_POST-PLACEBO_POST","PLACEBO_POST-IVIG_POST"),
                      group = c("IVIG_POST", "PLACEBO_POST"))
top_n_target = 250

#=============================================================
#Visualization of differential cell-cell interactions  
#ChordDiagram circos plots
#=============================================================
#get the top50 LR for circos plots
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 25, 
  rank_per_group = T
)
prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

#set the cluster colors
senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = col.cluster[senders_receivers]
colors_receiver =  col.cluster[senders_receivers]

#make the circos plots
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

jpeg(filename = paste0("circos plot IVIG_POST MNN top 50.jpeg"), width=1800 , height =1800,quality = 100,res = 300)
circos_list$IVIG_POST
dev.off()
svg(filename = paste0("circos plot IVIG_POST MNN top 50.svg"), width=10 , height =10)
circos_list$IVIG_POST
dev.off()

jpeg(filename = paste0("circos plot PLACEBO_POST MNN top 50.jpeg"), width=1800 , height =1800,quality = 100,res = 300)
circos_list$PLACEBO_POST
dev.off()
svg(filename = paste0("circos plot PLACEBO_POST MNN top 50.svg"), width=10 , height =10)
circos_list$PLACEBO_POST
dev.off()

jpeg(filename = paste0("circos plot legend.jpeg"), width=500 , height =1000,quality = 100,res = 300)
circos_list$legend
dev.off()
svg(filename = paste0("circos plot legend.svg"), width=4,height = 4)
circos_list$legend
dev.off()

#=============================================================
# quantify LR activity overall per celltype (top50)
#=============================================================  
# summarize changes in interaction between the conditions 
  MNN.df<-as.data.frame(multinichenet_output$prioritization_tables$group_prioritization_table_source)
  #filter for ligand or receptor siginficant change
  MNN.df<-MNN.df[MNN.df$p_val_ligand<0.05 | MNN.df$p_val_receptor<0.05,]
  #only take the IVIg direction (PLACEBO will simply be the opposit)
  MNN.df<-MNN.df[MNN.df$group=="IVIG_POST",]
  #reduce df
  MNN.df<-MNN.df[,c("group","sender","receiver","activity","id","direction_regulation")]
  #summerize
  MNN.df.new<-NULL
  for (direction in unique(MNN.df$direction_regulation)) {
    for (celltype in unique(MNN.df$sender)) {
      MNN.df.sub<-MNN.df[MNN.df$direction_regulation==direction,]
      interaction.count<-nrow(MNN.df.sub[MNN.df.sub$sender==celltype | MNN.df.sub$receiver==celltype,])
      activity_sum_S<-sum(MNN.df.sub$activity[MNN.df.sub$sender==celltype])
      activity_sum_R<-sum(MNN.df.sub$activity[MNN.df.sub$receiver==celltype])
      df<-data.frame("direction"=direction,"celltype"=celltype,"interaction.count"=interaction.count,"activity_sum_S"=activity_sum_S,"activity_sum_R"=activity_sum_R)
      MNN.df.new<-rbind(MNN.df.new,df)
    }}
  MNN.df.new$condition<-ifelse(MNN.df.new$direction=="up","IVIg","Placebo")
  MNN.df.new$condition<-factor( MNN.df.new$condition,c("Placebo","IVIg"))

  #bar plot grouped by comparison (all LR  with significant Ligand or receptor DE)
  p1<-ggplot(MNN.df.new, aes(x = celltype, y = activity_sum_R, fill = condition)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = col.ivig.condi[c(2,4)])+
    ylab("Receiver activity score")+xlab("")+ggtitle("Difference in LR interaction activity between IVIg and Placebo")+
    theme_classic()+
    theme(axis.text = element_text(size = 8),axis.text.x = element_blank(),axis.ticks.x = element_blank())
  
  p2<-ggplot(MNN.df.new, aes(x = celltype, y = activity_sum_S, fill = condition)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = col.ivig.condi[c(2,4)])+
    ylab("Sender activity score")+xlab("")+
    theme_classic()+
    theme(axis.text = element_text(size = 8),axis.text.x =element_text(angle = 45,hjust = 1))
  plot_grid(p1,p2,ncol = 1,align = "v")
  ggsave(filename = paste0("Sender and Receiver activity scores barplot grouped.jpeg"), width=3, height = 4)
  ggsave(filename = paste0("Sender and Receiver activity scores barplot grouped.svg"), width=3, height = 4)
  
# summarize changes in interaction between the conditions 
  #TOP 50
  prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, top_n = 25, rank_per_group = T)
  MNN.df<-as.data.frame(multinichenet_output$prioritization_tables$group_prioritization_table_source)
  #reduce to top 50
  MNN.df<-MNN.df[MNN.df$id %in% prioritized_tbl_oi_all$id,]
  #only take the IVIg direction (PLACEBO will simply be the opposit)
  MNN.df<-MNN.df[MNN.df$group=="IVIG_POST",]
  #reduce df
  MNN.df<-MNN.df[,c("group","sender","receiver","activity","id","direction_regulation")]
  #summerize
  MNN.df.new<-NULL
  for (direction in unique(MNN.df$direction_regulation)) {
    for (celltype in unique(MNN.df$sender)) {
      MNN.df.sub<-MNN.df[MNN.df$direction_regulation==direction,]
      interaction.count<-nrow(MNN.df.sub[MNN.df.sub$sender==celltype | MNN.df.sub$receiver==celltype,])
      activity_sum_S<-sum(MNN.df.sub$activity[MNN.df.sub$sender==celltype])
      activity_sum_R<-sum(MNN.df.sub$activity[MNN.df.sub$receiver==celltype])
      df<-data.frame("direction"=direction,"celltype"=celltype,"interaction.count"=interaction.count,"activity_sum_S"=activity_sum_S,"activity_sum_R"=activity_sum_R)
      MNN.df.new<-rbind(MNN.df.new,df)
    }}
  MNN.df.new$condition<-ifelse(MNN.df.new$direction=="up","IVIg","Placebo")
  MNN.df.new$condition<-factor( MNN.df.new$condition,c("Placebo","IVIg"))
  
  #bar plot grouped by comparison (all LR  with significant Ligand or receptor DE)
  p1<-ggplot(MNN.df.new, aes(x = celltype, y = activity_sum_R, fill = condition)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = col.ivig.condi[c(2,4)])+
    ylab("Receiver activity score")+xlab("")+ggtitle("Difference in LR interaction (Top50) activity between IVIg and Placebo")+
    theme_classic()+
    theme(axis.text = element_text(size = 8),axis.text.x = element_blank(),axis.ticks.x = element_blank())
  
  p2<-ggplot(MNN.df.new, aes(x = celltype, y = activity_sum_S, fill = condition)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = col.ivig.condi[c(2,4)])+
    ylab("Sender activity score")+xlab("")+
    theme_classic()+
    theme(axis.text = element_text(size = 8),axis.text.x =element_text(angle = 45,hjust = 1))
  plot_grid(p1,p2,ncol = 1,align = "v")
  ggsave(filename = paste0("Sender and Receiver activity scores top50 barplot grouped.jpeg"), width=3, height = 4)
  ggsave(filename = paste0("Sender and Receiver activity scores top50 barplot grouped.svg"), width=3, height = 4)
