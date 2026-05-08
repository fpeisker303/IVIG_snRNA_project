## Copyright (c) 2026 Fabian Peisker et al.

#IVIG project
#composition analysis

source("helper_functions_IVIG_publication_version.R")
source("IVIg2_snRNA_colors_publication.R")

#load final filtered and annotated snRNA object
sc_IVIg<-qread("IVIG_snRNA_harmony_f5.qs")
Idents(sc_IVIg)<-"celltypes.short"

#load function from speckle(propeller)
getTransformedProps <- function(clusters=clusters, sample=sample, transform=NULL){
  if(is.null(transform)) transform <- "logit"
  
  tab <- table(sample, clusters)
  props <- tab/rowSums(tab)
  if(transform=="asin"){
    message("Performing arcsin square root transformation of proportions")
    prop.trans <- asin(sqrt(props))
  }
  else if(transform=="logit"){
    message("Performing logit transformation of proportions")
    props.pseudo <- (tab+0.5)/rowSums(tab+0.5)
    prop.trans <- log(props.pseudo/(1-props.pseudo))
  }
  return(list(Counts=t(tab), TransformedProps=t(prop.trans), 
              Proportions=t(props)))
}

#pseudobulk function for concatenated meta data
sc_IVIg$three_groups_setup<-as.character(sc_IVIg$condition)
sc_IVIg$three_groups_setup[sc_IVIg$three_groups_setup=="PLACEBO_PRE"|sc_IVIg$three_groups_setup=="IVIG_PRE"]<-"PRE"
sc_IVIg$three_groups_setup<-factor(sc_IVIg$three_groups_setup,levels = c("PRE","PLACEBO_POST","IVIG_POST"))
ps<-create_pseudobulk_object(sc_IVIg,sample_ID = "orig.ident",cell_type_lable = "celltypes.short",add_ct_cs_faction = T)

#=========================================================================================
# major celltype lvl extended to more major ct (analysis using logit count)
# including adjustment for difference in cell coverage per sample !!!
#=========================================================================================
props.ct.ex<-getTransformedProps(clusters = sc_IVIg$celltypes.short.ext,sample = sc_IVIg$orig.ident)

#use limma to test for composition differences
group <- ps$condition
sample_cellnumber<- ps$sample_cellnumber
design <- model.matrix(~ 0 + group +sample_cellnumber)
fit1 <- lmFit(props.ct.ex$TransformedProps, design=design)

mycontr <- makeContrasts(groupIVIG_PRE - groupPLACEBO_PRE,
                         groupIVIG_POST - groupPLACEBO_POST,
                         groupIVIG_PRE - groupIVIG_POST,
                         groupPLACEBO_PRE - groupPLACEBO_POST,
                         diff_effect = (groupIVIG_PRE - groupIVIG_POST) - (groupPLACEBO_PRE - groupPLACEBO_POST),
                         levels=design)
fit2 <- contrasts.fit(fit1, mycontr)
fit2 <- eBayes(fit2)
topTable(fit2,coef = "diff_effect" , adjust="BH")
#test the relevant comparisons
res1<-topTable(fit2,coef = "groupIVIG_PRE - groupPLACEBO_PRE" , adjust="BH",number = Inf)
res1$clusters<-rownames(res1)
res1$group1<-"IVIG_PRE"
res1$group2<-"PLACEBO_PRE"
res1$contrast<-"groupIVIG_PRE - groupPLACEBO_PRE"
res2<-topTable(fit2,coef = "groupIVIG_POST - groupPLACEBO_POST" , adjust="BH",number = Inf)
res2$clusters<-rownames(res2)
res2$group1<-"IVIG_POST"
res2$group2<-"PLACEBO_POST"
res2$contrast<-"groupIVIG_POST - groupPLACEBO_POST"
res3<-topTable(fit2,coef = "groupIVIG_PRE - groupIVIG_POST" , adjust="BH",number = Inf)
res3$clusters<-rownames(res3)
res3$group1<-"IVIG_PRE"
res3$group2<-"IVIG_POST"
res3$contrast<-"groupIVIG_PRE - groupIVIG_POST"
res4<-topTable(fit2,coef = "groupPLACEBO_PRE - groupPLACEBO_POST" , adjust="BH",number = Inf)
res4$clusters<-rownames(res4)
res4$group1<-"PLACEBO_PRE"
res4$group2<-"PLACEBO_POST"
res4$contrast<-"groupPLACEBO_PRE - groupPLACEBO_POST"

#combine result to be added to the plot
res<-rbind(res1,res2,res3,res4)
rownames(res)<-NULL
res$clusters<-factor(res$clusters,levels= c("CM","EC","FIB","PC","SMC","Monocyte","Macrophage","cDC","MAST","Tcells","NK","NEURO","LYMEC","EndoEC","ADI"))
res$group1<-factor(res$group1,levels = c("PLACEBO_PRE","IVIG_PRE","IVIG_POST"))
res<-arrange(res,group1)
res<-arrange(res,clusters)

#add stat positions
#x positions
res$x <- as.numeric(factor(res$clusters, levels = unique(res$clusters)))
res$xmin <- res$x
res$xmax <- res$x
res$xmin[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]<-res$xmin[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]-0.3
res$xmax[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]<-res$xmax[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]+0.1
res$xmin[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]<-res$xmin[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]-0.1
res$xmax[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]<-res$xmax[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]+0.3
res$xmin[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]<-res$xmin[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]+0.1
res$xmax[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]<-res$xmax[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]+0.3
res$xmin[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]<-res$xmin[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]-0.3
res$xmax[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]<-res$xmax[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]-0.1
#y.position
y.position<-props.ct.ex$TransformedProps %>%
  melt() %>%
  group_by(clusters) %>%
  summarize(y.position = max(value, na.rm = TRUE))
res<-left_join(res,y.position,by="clusters")
values_to_add <- c(0.1, 0.2, 0.3,0.4)
rows_to_modify <- seq(2, nrow(res), by = 4)  # Identify starting rows for each block
for (i in 1:length(values_to_add)) {
  current_rows <- rows_to_modify + (i - 1)
  current_rows <- current_rows[current_rows <= nrow(res)]  # Ensure no out-of-bound rows
  res$y.position[current_rows] <- res$y.position[current_rows] + values_to_add[i]
}
#sig style
res<-arrange(res,adj.P.Val)
res$p.adj.signif.NA<-res$adj.P.Val
res$p.adj.signif.NA[res$adj.P.Val > 0.1]<-NA
res$p.adj.signif.NA[res$adj.P.Val < 0.1]<-'#'
res$p.adj.signif.NA[res$adj.P.Val < 0.05]<-'*'
res$p.adj.signif.NA[res$adj.P.Val < 0.01]<-'**'

#to plot the transformed proportions
breakdown.df<-melt(props.ct.ex$TransformedProps)
sample <- ps$orig_ident
condition<-ps$condition
df<-data.frame(sample,condition)
breakdown.df<-left_join(breakdown.df,df,by = "sample")
breakdown.df$clusters<-factor(breakdown.df$clusters,levels = c("CM","EC","FIB","PC","SMC","Monocyte","Macrophage","cDC","MAST","Tcells","NK","NEURO","LYMEC","EndoEC","ADI"))

# Unique data frame to draw one rectangle per facet only
facet_backgrounds<-data.frame(x=1:15,clusters=levels(res$clusters))
facet_backgrounds$xmin<-facet_backgrounds$x-0.5
facet_backgrounds$xmax<-facet_backgrounds$x+0.5
facet_backgrounds$ymin<--Inf
facet_backgrounds$ymax<-Inf
facet_backgrounds$clusters<-factor(facet_backgrounds$clusters,levels = c("CM","EC","FIB","PC","SMC","Monocyte","Macrophage","cDC","MAST","Tcells","NK","NEURO","LYMEC","EndoEC","ADI"))

# Define custom colors for each facet
facet_colors <- c("CM" = "#88CCEE4D", "EC" = "#DDCC774D", "FIB" = "#CC66774D",
                  "PC"="#1177334D","SMC"="#1177334D",
                  "Monocyte"="#3322884D","Macrophage"="#3322884D","cDC"="#3322884D","MAST"="#3322884D",
                  "Tcells"="#AA44994D","NK"="#AA44994D",
                  "EndoEC"="#44AA994D","LYMEC"="#9999334D","NEURO"="#8888884D","ADI"="#6611004D")

#no points, bar plot version
#shift to positive scale
shift_by<-abs(min(breakdown.df$value))
breakdown.df$value<-breakdown.df$value + shift_by
#Calculate mean and SD per group
summary_df <- breakdown.df %>%
  group_by(clusters,condition) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )
#adjust stats position
res$y.position<-res$y.position + shift_by

ggplot(summary_df,  aes(x = clusters, y = mean)) +
  geom_rect(data = facet_backgrounds, aes(fill = clusters, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE) +
  scale_fill_manual(values = gsub("4D","20",facet_colors))+
  ggnewscale::new_scale_fill()+
  geom_bar(stat = "identity",aes(fill= condition),position=position_dodge2(padding=0.1))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position = position_dodge2(),linewidth=0.2) +
  scale_fill_manual(values =col.ivig.condi)+
  ylab("Relative Celltype Abundance")+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15))+
  scale_y_continuous(expand = c(0, 0)) +
  stat_pvalue_manual(res, label = "p.adj.signif.NA", tip.length = 0.05,fontface="bold",size=4)

ggsave(filename = "logit ct_ext count per condition with limma stat_BAR.jpeg", width=10 , height = 6)
ggsave(filename = "logit ct_ext count per condition with limma stat_BAR.svg", width=10 , height = 6)

#=========================================================================================
# testing only the CM subtypes
# normalization to total nuc count per sample
# including adjustment for difference in cell coverage per sample !!!
#=========================================================================================

#load
sc_IVIg.CM<-qread("subclustering of CM_annotated.qs")

#to get total get ps on full data
total_count<-ps$sample_cellnumber

#get myeloid numbers
tab <- table(sc_IVIg.CM$subcelltypes.short, sc_IVIg.CM$orig.ident)
total_count<-total_count[colnames(tab)]

#add pseudocount
tab<-tab+0.5
total_count<-total_count+(nrow(tab)*0.5)
# Divide each column by the corresponding vector value
props.pseudo <- sweep(tab, 2, total_count, FUN = "/")
#logit
prop.CM <- log(props.pseudo/(1-props.pseudo))

#use limma to test for composition differences
group <- ps$condition
sample_cellnumber<- ps$sample_cellnumber
design <- model.matrix(~ 0 + group +sample_cellnumber)
fit1 <- lmFit(prop.CM, design=design)

mycontr <- makeContrasts(groupIVIG_PRE - groupPLACEBO_PRE,
                         groupIVIG_POST - groupPLACEBO_POST,
                         groupIVIG_PRE - groupIVIG_POST,
                         groupPLACEBO_PRE - groupPLACEBO_POST,
                         diff_effect = (groupIVIG_PRE - groupIVIG_POST) - (groupPLACEBO_PRE - groupPLACEBO_POST),
                         levels=design)
fit2 <- contrasts.fit(fit1, mycontr)
fit2 <- eBayes(fit2)
topTable(fit2,coef = "diff_effect" , adjust="BH")
#test the relevant comparisons
res1<-topTable(fit2,coef = "groupIVIG_PRE - groupPLACEBO_PRE" , adjust="BH",number = Inf)
res1$clusters<-rownames(res1)
res1$group1<-"IVIG_PRE"
res1$group2<-"PLACEBO_PRE"
res1$contrast<-"groupIVIG_PRE - groupPLACEBO_PRE"
res2<-topTable(fit2,coef = "groupIVIG_POST - groupPLACEBO_POST" , adjust="BH",number = Inf)
res2$clusters<-rownames(res2)
res2$group1<-"IVIG_POST"
res2$group2<-"PLACEBO_POST"
res2$contrast<-"groupIVIG_POST - groupPLACEBO_POST"
res3<-topTable(fit2,coef = "groupIVIG_PRE - groupIVIG_POST" , adjust="BH",number = Inf)
res3$clusters<-rownames(res3)
res3$group1<-"IVIG_PRE"
res3$group2<-"IVIG_POST"
res3$contrast<-"groupIVIG_PRE - groupIVIG_POST"
res4<-topTable(fit2,coef = "groupPLACEBO_PRE - groupPLACEBO_POST" , adjust="BH",number = Inf)
res4$clusters<-rownames(res4)
res4$group1<-"PLACEBO_PRE"
res4$group2<-"PLACEBO_POST"
res4$contrast<-"groupPLACEBO_PRE - groupPLACEBO_POST"

#combine result to be added to the plot
res<-rbind(res1,res2,res3,res4)
rownames(res)<-NULL
res$clusters<-factor(res$clusters,levels= levels(sc_IVIg.CM$subcelltypes.short))
res$group1<-factor(res$group1,levels = c("PLACEBO_PRE","IVIG_PRE","IVIG_POST"))
res<-arrange(res,group1)
res<-arrange(res,clusters)

#add stat positions
#x positions
res$x <- as.numeric(factor(res$clusters, levels = unique(res$clusters)))
res$xmin <- res$x
res$xmax <- res$x
res$xmin[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]<-res$xmin[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]-0.3
res$xmax[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]<-res$xmax[res$contrast=="groupIVIG_PRE - groupPLACEBO_PRE"]+0.1
res$xmin[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]<-res$xmin[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]-0.1
res$xmax[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]<-res$xmax[res$contrast=="groupIVIG_POST - groupPLACEBO_POST"]+0.3
res$xmin[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]<-res$xmin[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]+0.1
res$xmax[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]<-res$xmax[res$contrast=="groupIVIG_PRE - groupIVIG_POST"]+0.3
res$xmin[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]<-res$xmin[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]-0.3
res$xmax[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]<-res$xmax[res$contrast=="groupPLACEBO_PRE - groupPLACEBO_POST"]-0.1
#y.position
y.position<-prop.CM %>%
  melt() %>%
  group_by(Var1) %>%
  summarize(y.position = max(value, na.rm = TRUE))
colnames(y.position)[1]<-"clusters"
res<-left_join(res,y.position,by="clusters")
values_to_add <- c(0.1, 0.2, 0.3,0.4)
rows_to_modify <- seq(2, nrow(res), by = 4)  # Identify starting rows for each block
for (i in 1:length(values_to_add)) {
  current_rows <- rows_to_modify + (i - 1)
  current_rows <- current_rows[current_rows <= nrow(res)]  # Ensure no out-of-bound rows
  res$y.position[current_rows] <- res$y.position[current_rows] + values_to_add[i]
}
#sig style
res<-arrange(res,P.Value)
res$p.adj.signif.NA<-res$P.Value
res$p.adj.signif.NA[res$P.Value > 0.1]<-NA
res$p.adj.signif.NA[res$P.Value < 0.1]<-'#'
res$p.adj.signif.NA[res$P.Value < 0.05]<-'*'
res$p.adj.signif.NA[res$P.Value < 0.01]<-'**'

#to plot the transformed proportions
breakdown.df<-melt(prop.CM)
colnames(breakdown.df)<-c("clusters","sample","value")
sample <- ps$orig_ident
condition<-ps$condition
df<-data.frame(sample,condition)
breakdown.df<-left_join(breakdown.df,df,by = "sample")
breakdown.df$clusters<-factor(breakdown.df$clusters,levels = levels(sc_IVIg.CM$subcelltypes.short))

# Unique data frame to draw one rectangle per facet only
facet_backgrounds<-data.frame(x=1:length(levels(res$clusters)),clusters=levels(res$clusters))
facet_backgrounds$xmin<-facet_backgrounds$x-0.5
facet_backgrounds$xmax<-facet_backgrounds$x+0.5
facet_backgrounds$ymin<--Inf
facet_backgrounds$ymax<-Inf
facet_backgrounds$clusters<-factor(facet_backgrounds$clusters,levels = levels(sc_IVIg.CM$subcelltypes.short))

# Define custom colors for each facet
facet_colors <- paste0(col.CM,"4D")
names(facet_colors)<-levels(sc_IVIg.CM$subcelltypes.short)

ggplot(breakdown.df, aes(x=clusters, y=value)) +
  geom_rect(data = facet_backgrounds, aes(fill = clusters, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE) +
  scale_fill_manual(values = facet_colors)+
  ggnewscale::new_scale_fill()+
  geom_boxplot(aes(fill= condition),position=position_dodge2(padding=0.1), pch=21, show.legend = T,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.8),pch=21,size=1, aes(fill= condition), show.legend = T)+
  scale_fill_manual(values =col.ivig.condi)+
  ylab("Relative Celltype Proportions [logit(rel. count)]")+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15))+
  stat_pvalue_manual(res, label = "p.adj.signif.NA", tip.length = 0.01,fontface="bold",size=4)
ggsave(filename = "logit CM cs count per condition with limma stat_norm total count.jpeg", width=8 , height = 6)
ggsave(filename = "logit CM cs count per condition with limma stat_norm total count.svg", width=8 , height = 6)


