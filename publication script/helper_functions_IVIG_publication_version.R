## Copyright (c) 2026 Fabian Peisker et al.

#helper functions used in the IVIG project
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(clustree)
library(writexl)
library(dplyr)
library(ggpubr)
library(scales)
library(dittoSeq)
library(reshape2)
library(readxl)
library(stringr)
library(scDblFinder)
library(SoupX)
library(harmony)
library(DESeq2)
library(Matrix.utils)
library(pheatmap)
library(biomaRt)
library(msigdbr)
library(SCP)
library(GSVA)
library(limma)
library(qs)
library(edgeR)

#set seurat to work with v3 objects
options(Seurat.object.assay.version = "v3")

#---fix seed
set.seed("123")

#--------SoupX function------------------------------
run.SoupX <- function(sc.object,path_cellranger_out){
  
  # Load data and estimate soup profile
  toc = Read10X(paste0(path_cellranger_out,"/filtered_feature_bc_matrix/"))
  clusters <- sc.object$seurat_clusters
  names(clusters) <- gsub(paste0(unique(sc.object$orig.ident),"_"),"",names(clusters))
  names(clusters) <-gsub("_","",names(clusters))
  #reduce filtered input to whats in the sc.object
  toc<-toc[,colnames(toc) %in% names(clusters)]
  
  tod = Read10X(paste0(path_cellranger_out,"/raw_feature_bc_matrix/"))
  sc = SoupChannel(tod, toc)
  
  #add clusters from seurat and UMAP
  sc = setClusters(sc,clusters)
  seurat.umap<-Embeddings(sc.object,reduction = "umap")
  row.names(seurat.umap)<-gsub(paste0(unique(sc.object$orig.ident),"_"),"",row.names(seurat.umap))
  row.names(seurat.umap)<-gsub("_","",row.names(seurat.umap))
  sc = setDR(sc, seurat.umap)
  
  #calculate correction
  sc <- autoEstCont(sc,forceAccept = T,doPlot = F)
  
  #adjust count matrix
  out <- adjustCounts(sc)
  
  #generate new seurat object based on adjusted count matrix
  out.sc <- CreateSeuratObject(counts=out, project=sc.object@project.name,min.cells = 3, min.features = 200)
  out.sc <- NormalizeData(out.sc, verbose = FALSE)
  out.sc <- FindVariableFeatures(out.sc, verbose = FALSE)
  out.sc <- recluster_RNA(out.sc,0.5)
  
  #add UMAP from the non-adjusted for comparison
  old.umap<-Embeddings(sc.object,reduction = "umap")
  row.names(old.umap)<-gsub(paste0(unique(sc.object$orig.ident),"_"),"",row.names(old.umap))
  row.names(old.umap)<-gsub("_","",row.names(old.umap))
  old.umap<-old.umap[row.names(old.umap) %in% colnames(out.sc),]
  out.sc[["oldumap"]] <- CreateDimReducObject(embeddings = old.umap,key = "oldumap_",assay = "RNA")
  
  #plot the top5 soup genes in old and new umap
  #get top5 genes that are probably soup
  top20.soup <- rownames(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20))
  top5.in.the.data <-top20.soup[top20.soup %in% row.names(out.sc)][1:5]
  pg.list=NULL
  for (top5 in top5.in.the.data) {
    f1<-plotChangeMap(sc, out, top5)+theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "none",title = element_text(size = 8))
    f2<-FeaturePlot(sc.object,features = top5)+theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "none",title = element_text(size = 8))
    f3<-FeaturePlot(out.sc,features = top5,reduction = "oldumap")+theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "none",title = element_text(size = 8))
    f4<-FeaturePlot(out.sc,features = top5,reduction = "umap")+theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "none",title = element_text(size = 8))
    pg.list[[top5]]<-plot_grid(plotlist = list(f1,f2,f3,f4),nrow = 1)
    
  }
  dev.off()
  
  
  return(list(out.sc,pg.list))#returns the sc object based on the filtered matrix and the plots of the top5
}


#------function for quality check of unfiltered scRNA dataset-----
assess_quality <- function(sc.object,res=0.5,runSoupX=T,sample_path,upperFeatureCutoff = 4000,lowerFeatureCutoff=500,mtCutoff =10){ 
  #partitial source https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html
  #sc.object is the seurat object
  #res is the resolution used for clustering. default 0.5
  #sample_path is the path to the cellranger out. required only if runSoupX = true
  
  #step1 normalize and perform standart seurat clustering pipeline, SoupX and scDblFinder require clusters
  sc.object <- NormalizeData(sc.object, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  sc.object <- FindVariableFeatures(sc.object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sc.object <- recluster_RNA(sc.object,res)
  
  #step2 if required, run SoupX to correct for ambientRNA (deflaut settings of SoupX)
  if (runSoupX) {
    require(SoupX)
    SoupX.result <- run.SoupX(sc.object,path_cellranger_out = sample_path)
    #continue with the SoupX corrected matrix
    sc.object <-SoupX.result[[1]]
  } 
  DefaultAssay(sc.object)="RNA"
  
  #step3 check for mitrochondial and ribosomal genes
  sc.object[['percent.mt']] = PercentageFeatureSet(sc.object, pattern = '^MT-')
  sc.object[['percent.rb']] = PercentageFeatureSet(sc.object, pattern = '^RP[SL][[:digit:]]')
  
  #step4 score for doublets
  doublets <- scDblFinder(sce = as.SingleCellExperiment(sc.object), clusters = sc.object$seurat_clusters)
  sc.object$doublet_score <- doublets$scDblFinder.score
  sc.object$doublet <- doublets$scDblFinder.class
  
  #step5 add quality assessment to metadata
  #doublet?
  sc.object[['QC']] <- ifelse(sc.object$doublet == 'singlet','Pass','Doublet')
  #low nFeature?
  sc.object[['QC']] <- ifelse(sc.object$nFeature_RNA < lowerFeatureCutoff,'Low_nFeature',sc.object$QC)
  #high nFeature?
  sc.object[['QC']] <- ifelse(sc.object$nFeature_RNA > upperFeatureCutoff,"High_nFeature",sc.object$QC)
  #high mt% ?
  sc.object[['QC']] <- ifelse(sc.object$percent.mt > mtCutoff,'High_MT',sc.object$QC)
  #table(sc.object[['QC']])
  
  #step6 plot overview of the quality and returned (unfiltered) seurat object
  v1 <- VlnPlot(sc.object,group.by = "orig.ident",features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt','percent.rb'), ncol = 4, pt.size = 0,combine = F)
  v1<-lapply(v1,function(x){x<-x+NoLegend()+theme(text = element_text(size = 4),axis.text.y = element_text(size = 4),axis.title.x = element_blank(),axis.text.x = element_blank())})
  v1<-plot_grid(plotlist = v1,nrow = 1)
  f1 <- FeaturePlotCombined(sc.object,features = c('nFeature_RNA'),pt.size = 0.2,order = T,title.size = 8,nolegend = T)
  #f2 <- FeaturePlotCombined(sc.object,features = c('nCount_RNA'),pt.size = 0.2,order = T,title.size = 8) 
  f3 <- FeaturePlotCombined(sc.object,features = c('percent.mt'),pt.size = 0.2,order = T,title.size = 8,nolegend = T) 
  f4 <- FeaturePlotCombined(sc.object,features = c('percent.rb'),pt.size = 0.2,order = T,title.size = 8,nolegend = T)
  d1 <- DimPlot(sc.object,group.by = "QC",pt.size = 0.2)+theme(title = element_blank(),axis.text = element_blank(),text=element_text(size=8),
                                                               legend.key.height= unit(0.1, 'cm'),
                                                               legend.key.width= unit(0.1, 'cm'),
                                                               legend.position = c(0.6,0.9))
  #put togehter
  pl1 <- ggarrange(v1,f1,f3,f4,nrow = 2,ncol=2)
  pl2 <- ggarrange(pl1,d1,nrow = 1)
  #add soupx results
  if(runSoupX){
    pl3 <- SoupX.result[[2]]
    pg2<-plot_grid(plotlist = list(pl3[[1]],pl3[[2]],pl3[[3]],pl3[[4]],pl3[[5]],pl2),ncol = 1,rel_heights = c(1,1,1,1,1,2))
  }else{
    pg2<-pl2
  }
  
  
  
  return(list(sc.object,pg2))
}


#----function to perfrom clustering following the default seurat pipeline------
recluster_RNA <- function(object,res){
  DefaultAssay(object)="RNA"
  # Run the standard workflow for visualization and clustering
  object <- NormalizeData(object,verbose=F)
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:30,repulsion.strength = 5, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  print(DimPlot(object))
  message("plot not saved")
  return(object)
}


#-------optimized nice function to plot combination of feature plots--------------
FeaturePlotCombined<-function(sc.object,features,pt.size=1,order=T, saveIndividual=F,title.size=30,reduction="umap",clusters="seurat_clusters",returnPLot=F,addClustPlot=F, nolegend=F, slot="scale.data",raster=F,assay="RNA"){
  gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
  featurePlotList=NULL
  DefaultAssay(sc.object)=assay
  if(saveIndividual){dir.create("individual_feature_plots")}
  if(nolegend){
    for (feature in features) {
      if (feature %in% rownames(sc.object@assays$RNA)){
        sc.object<-ScaleData(sc.object,features=feature,verbose=F)
        featurePlotList[[feature]]=FeaturePlot(sc.object, features = feature, min.cutoff = "q9",pt.size = pt.size,order = order,reduction = reduction, slot=slot,raster = raster) + 
          scale_colour_gradientn(colours =gradient.col)+guides(colour="none")+
          theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))
        if (reduction!="umap"){featurePlotList[[feature]]=featurePlotList[[feature]]+xlim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,1]),max(sc.object@reductions[[reduction]]@cell.embeddings[,1])))+
          ylim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,2]),max(sc.object@reductions[[reduction]]@cell.embeddings[,2])))}
        if(saveIndividual){print(featurePlotList[[feature]]);ggsave(filename = paste0("individual_feature_plots/",feature,".jpeg"), width=10 , height = 10)}
      }else{message(paste0("feature ",feature," not in the gene data!"))}
      if (feature %in% colnames(sc.object[[]])) {
        featurePlotList[[feature]]=FeaturePlot(sc.object, features = feature,pt.size = pt.size,order = order,reduction = reduction, slot=slot,raster = raster) + 
          scale_colour_gradientn(colours =gradient.col)+guides(colour="none")+
          theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))
        if (reduction!="umap"){featurePlotList[[feature]]=featurePlotList[[feature]]+xlim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,1]),max(sc.object@reductions[[reduction]]@cell.embeddings[,1])))+
          ylim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,2]),max(sc.object@reductions[[reduction]]@cell.embeddings[,2])))}
        if(saveIndividual){print(featurePlotList[[feature]]);ggsave(filename = paste0("individual_feature_plots/",feature,".jpeg"), width=10 , height = 10)}
      }
    }#end of for
    if (addClustPlot) {featurePlotList[["clusters"]]=DimPlot(sc.object, group.by = clusters,reduction = reduction, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))}
    if (returnPLot){return(ggarrange(plotlist = featurePlotList,common.legend = T,legend = "right"))}else{
      print(ggarrange(plotlist = featurePlotList,common.legend = T,legend = "right"))}
    
  }else{
    for (feature in features) {
      if (feature %in% rownames(sc.object@assays$RNA)){
        sc.object<-ScaleData(sc.object,features=feature,verbose=F)
        featurePlotList[[feature]]=FeaturePlot(sc.object, features = feature, min.cutoff = "q9",pt.size = pt.size,order = order,reduction = reduction, slot=slot,raster = raster) + 
          scale_colour_gradientn(colours =gradient.col)+
          theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))
        if (reduction!="umap"){featurePlotList[[feature]]=featurePlotList[[feature]]+xlim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,1]),max(sc.object@reductions[[reduction]]@cell.embeddings[,1])))+
          ylim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,2]),max(sc.object@reductions[[reduction]]@cell.embeddings[,2])))}
        if(saveIndividual){print(featurePlotList[[feature]]);ggsave(filename = paste0("individual_feature_plots/",feature,".jpeg"), width=10 , height = 10)}
      }else{message(paste0("feature ",feature," not in the gene data!"))}
      if (feature %in% colnames(sc.object[[]])) {
        featurePlotList[[feature]]=FeaturePlot(sc.object, features = feature,pt.size = pt.size,order = order,reduction = reduction, slot=slot,raster = raster) + 
          scale_colour_gradientn(colours =gradient.col)+
          theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))
        if (reduction!="umap"){featurePlotList[[feature]]=featurePlotList[[feature]]+xlim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,1]),max(sc.object@reductions[[reduction]]@cell.embeddings[,1])))+
          ylim(c(min(sc.object@reductions[[reduction]]@cell.embeddings[,2]),max(sc.object@reductions[[reduction]]@cell.embeddings[,2])))}
        if(saveIndividual){print(featurePlotList[[feature]]);ggsave(filename = paste0("individual_feature_plots/",feature,".jpeg"), width=10 , height = 10)}
      }
    }#end of for
    if (addClustPlot) {featurePlotList[["clusters"]]=DimPlot(sc.object, group.by = clusters,reduction = reduction, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()+theme(axis.title = element_blank(),axis.text = element_blank(),plot.title = element_text(size=title.size))}
    if (returnPLot){return(ggarrange(plotlist = featurePlotList,common.legend = T,legend = "right"))}else{
      print(ggarrange(plotlist = featurePlotList,common.legend = T,legend = "right"))}
  }#end of no legend
}#end of function

#---------function to generate a seurat object with metacell/pseudobulks instead of cells------------
#---------function to generate a seurat object with metacell/pseudobulks instead of cells------------
create_pseudobulk_object<-function(seurat_object,sample_ID,return_seurat=T,cell_type_lable=NULL,cell_states_label=NULL,
                                   average_meta_feature=NULL,pseudobulk_per_ct=F,add_ct_cs_faction=F,cell_number_cutoff=20){
  #set group labels for the pseudobulking, either by sample_ID, or by sample_ID + ct label to have pseudobulks per ct
  if (pseudobulk_per_ct) {
    if (is.null(cell_type_lable)) { stop("for pseudobulking separated per celltype, cell_type_lable parameter must be set")}
    seurat_object$sample_ID_ct<-paste0(seurat_object[[cell_type_lable]][,1],"_",  seurat_object[[sample_ID]][,1])
    groups<-FetchData(seurat_object,vars = "sample_ID_ct")
    sample_ID<-"sample_ID_ct"
  }else{
    groups<-FetchData(seurat_object,vars = sample_ID)
  }
  
  #exclude samples with less than 20 cells contributing to a pseudobulk
  groups_keep<-names(table(groups))[table(groups)>cell_number_cutoff]
  groups_removed<-names(table(groups))[table(groups)<cell_number_cutoff]
  if (length(groups_removed)!=0) {
    Idents(seurat_object)<-sample_ID
    seurat_object<-subset(seurat_object,idents=groups_keep)
    message(paste0("Samples ",groups_removed," have less then ",cell_number_cutoff," cells contributing and are therefore excluded from pseudobulking \n"))
  } else {
    message("All samples have a sufficient number of cell to generate a pseudobulk \n")
  }
  
  #generate pseudobulks and create a metacell seurat object
  groups<-FetchData(seurat_object,vars = sample_ID)
  message(paste0("generating ",length(unique(groups[,1]))," metacells/pseudobulks"))
  average.exp <- aggregate.Matrix(t(seurat_object@assays$RNA@counts), 
                                  groupings = groups[,sample_ID], fun = "sum")
  #check if values are not integer and report (counts should usually be integers)
  message("First 100 gene expression values of the pseudobulk count matrix, before rounding eventually float values")
  # Sample matrix with both integers and floats
  mat <- average.exp[1:nrow(average.exp),1:100]
  print(mat)
  
  #get ride of genes with only 0 in all groups
  average.exp<-round(t(average.exp))
  average.exp<-average.exp[rowSums(average.exp[])>0,]
  
  #create metadata
  #automatically get all meta where the individual sample has only one unique value
  #other meta data cant be considered
  meta.data <- FetchData(seurat_object,vars = colnames(seurat_object@meta.data))
  #save table with cellnumber per sample
  sample.cellnumber<-table(meta.data[,sample_ID])
  sample.cellnumber<-data.frame("cell_number"=c(sample.cellnumber))
  sample.cellnumber<-sample.cellnumber[colnames(average.exp),]#rm excluded samples
  #rm all meta.data with more values than sample number and reduce to a small meta.data df
  df<-NULL
  for (sample in unique(meta.data[,sample_ID])) {
    df[[sample]]<-sapply(meta.data[meta.data[,sample_ID]==sample,], function(x) length(unique(x))>1 ) }
  meta.data<-meta.data[,rowSums(as.data.frame(df))==0]
  meta.data<-meta.data %>% distinct()
  rownames(meta.data)<-meta.data[,sample_ID]
  #rm excluded samples
  meta.data<-meta.data[colnames(average.exp),]
  #combine
  meta.data<-cbind(meta.data,sample.cellnumber)
  
  #add cell type numbers if annotation is added
  if (!is.null(cell_type_lable)) {
    message(paste0("adding per sample cell type number"))
    cell_type.number<-FetchData(seurat_object,vars=c(sample_ID,cell_type_lable))
    cell_type.number<-as.data.frame(table(cell_type.number))
    cell_type.number<-dcast(as.data.frame(cell_type.number),as.formula(paste0(colnames(cell_type.number)[1], "~" ,colnames(cell_type.number)[2])),value.var = "Freq")
    rownames(cell_type.number)<-cell_type.number[,1]
    cell_type.number<-cell_type.number[,-1]
    colnames(cell_type.number)<-gsub(" ","_",colnames(cell_type.number))
    colnames(cell_type.number)<-paste0("ct_",colnames(cell_type.number)) #add common identifier label of cell type number
    #only overlapping samples
    cell_type.number<-cell_type.number[rownames(meta.data),]
    #combine
    meta.data<-cbind(meta.data,cell_type.number)
    #optional: also add the faction of the cell types per sample (sometime more informative)
    if (add_ct_cs_faction) {
      #use dittobarplot to easy get the fractions
      df_ct_pct<-dittoBarPlot(seurat_object,var = cell_type_lable,group.by = sample_ID)
      df_ct_pct<-dcast(df_ct_pct$data, grouping ~ label,value.var = "percent")
      row.names(df_ct_pct)<-df_ct_pct$grouping
      df_ct_pct$grouping<-NULL
      colnames(df_ct_pct)<-paste0("ct_",colnames(df_ct_pct),"_pct")
      df_ct_pct<-df_ct_pct[rownames(meta.data),] #remove data this is not needed and align order
      meta.data<-cbind(meta.data,df_ct_pct)
    }
  }
  
  #add cell status number if annotation is added
  if (!is.null(cell_states_label)) {
    message(paste0("adding per sample cell state number"))
    cell_state.number<-FetchData(seurat_object,vars=c(sample_ID,cell_states_label))
    cell_state.number<-as.data.frame(table(cell_state.number))
    cell_state.number<-dcast(as.data.frame(cell_state.number),as.formula(paste0(colnames(cell_state.number)[1], "~" ,colnames(cell_state.number)[2])),value.var = "Freq")
    rownames(cell_state.number)<-cell_state.number[,1]
    cell_state.number<-cell_state.number[,-1]
    colnames(cell_state.number)<-gsub(" ","_",colnames(cell_state.number))
    colnames(cell_state.number)<-paste0("cs_",colnames(cell_state.number)) #add common identifier label of cell type number
    #only overlapping samples
    cell_state.number<-cell_state.number[rownames(meta.data),]
    #combine
    meta.data<-cbind(meta.data,cell_state.number)
    #optional: also add the faction of the cell state per sample (sometime more informative)
    if (add_ct_cs_faction) {
      #for cell state factions normalize to number of high level celltype
      
      #Endothelial cell types make some problems due to some EndoEC and LymEC clustering different between ct and cs
      #manual correction in case those ct and cs are included
      annotations<-unique(c(as.character(seurat_object[[cell_states_label]][,1]),as.character(seurat_object[[cell_type_lable]][,1])))
      if(all(c("EndoEC","LymEC","LYMEC") %in% annotations)){
        message("problematic cell type and state label for EndoEC and LymEC found")
        #get ct and cs relation
        df_cell_anno <- seurat_object@meta.data[c(cell_states_label,cell_type_lable)]
        df_cell_anno <- split(as.character(df_cell_anno$subcelltypes.short), df_cell_anno$celltypes.short)
        df_cell_anno <- lapply(df_cell_anno, unique)
        #correct for EC
        df_cell_anno$EndoEC<-c("EndoEC")
        df_cell_anno$LYMEC<-c("LymEC")
        df_cell_anno$EC<-c("capillaryEC", "angioEC",  "veinEC" , "arteryEC" )
      }else{
        df_cell_anno<-df_cell_anno[unique(melt(df_cell_anno)[,"L1"])]
      }
      #isolate relevant number and calculate freq
      breakdown<-table(seurat_object[[c(cell_states_label,sample_ID)]])
      for (ct in names(df_cell_anno)) {
        if (ct==names(df_cell_anno)[1]) {
          df <- breakdown[rownames(breakdown) %in% df_cell_anno[[ct]],]
          combined_df<-apply(df, 2, function(x){(x/sum(x))*100 })
        }else{
          df <- breakdown[rownames(breakdown) %in% df_cell_anno[[ct]],]
          if (is.null(nrow(df))) {
            message(paste0("cell state: ",ct," is the only state of that type, therefore normalized to total cell number"))
            df<-t(as.data.frame((df / colSums(breakdown))*100))
            rownames(df) <- df_cell_anno[[ct]]
          } else {
            df<-apply(df, 2, function(x){(x/sum(x))*100 })
          }
          combined_df<-rbind(combined_df,df)
        }
      }#end of for df_cell_anno
      #add to meta data
      df_cs_pct<-t(combined_df)
      colnames(df_cs_pct)<-paste0("cs_",colnames(df_cs_pct),"_pct")
      df_cs_pct<-df_cs_pct[rownames(meta.data),] #remove data this is not needed and align order
      meta.data<-cbind(meta.data,df_cs_pct)
    }
  }
  
  #for continues meta features an average per pseudobulk can be calculated and added to meta data
  if (!is.null(average_meta_feature)) {
    message(paste0(c("adding per sample average for: ",average_meta_feature)))
    number.data<-FetchData(seurat_object,vars = c(sample_ID,average_meta_feature))
    colnames(number.data)[1]<-"sample_ID"
    averages <- number.data%>%
      group_by(sample_ID) %>%
      summarize_all(mean) %>% as.data.frame()
    rownames(averages)<-averages$sample_ID
    averages$sample_ID<-NULL
    colnames(averages)<-paste0("avg_",colnames(averages))
    #only overlapping samples
    averages<-averages[rownames(meta.data),]
    meta.data<-cbind(meta.data,averages)
  }
  
  #adjust problematic colnames of metadata
  colnames(meta.data)<-gsub("-","_",colnames(meta.data))
  colnames(meta.data)<-gsub("[+]","pos",colnames(meta.data))
  colnames(meta.data)<-gsub("[.]","_",colnames(meta.data))
  colnames(meta.data)<-gsub("[/]","_",colnames(meta.data))
  
  #remove factor levels that dont exist anymore
  meta.data <- meta.data %>%
    mutate(across(where(is.factor), ~ factor(.)))
  
  if (return_seurat) {
    #create meta seurat object
    options(Seurat.object.assay.version = "v3")
    seurat_object.meta <- CreateSeuratObject(counts = average.exp,meta.data = meta.data)
    seurat_object.meta@meta.data$orig.ident<-NULL
    return(seurat_object.meta)
  }else{
    #return the pseudobulk count matrix and meta.data dataframe as list
    message("return the pseudobulk count matrix and meta.data dataframe as list")
    return(list("count_matrix"=average.exp,"meta_data"=meta.data))
  }
  
}


#---------filter out genes that are very low expressed and only in a very few samples covered-----------
pseudobulk.QC.filter<-function(pseudo.object,min.sample=NULL,min.percent.sample=NULL,use_edgeR=F){
  if (!is.null(min.percent.sample)&!is.null(min.sample)) { #this check does not work properly!
    message("filtering can only be either by min.sample cut off or min.percent.sample. not both.")
    stop()}
  if (use_edgeR) {
    dge <- DGEList(counts = as.matrix(pseudo.object@assays$RNA@counts))
    keep <- filterByExpr(dge)
    genes_to_keep<-names(keep)[keep]
  }else{
    #number of samples vs log10 total sum of counts
    how_often_in_sample <- rowSums(!pseudo.object@assays$RNA$counts == 0)
    names(how_often_in_sample) <- rownames(pseudo.object@assays$RNA$counts)
    if (!is.null(min.sample)) {
      #genes below threshhold of min.sample
      genes_to_keep<-names(how_often_in_sample)[!how_often_in_sample<min.sample]
    }
    if (!is.null(min.percent.sample)) {
      #check how often genes are expressed in samples overall
      #remove genes that are expressed in less than min.percent.sample samples.
      min.sample<-round(max(how_often_in_sample)*min.percent.sample)
      message(paste0("using min.percent.sample: ",min.percent.sample," gives a cut off value of: ",min.sample))
      genes_to_keep<-names(how_often_in_sample)[!how_often_in_sample<min.sample]
    }
  }
  #display how many genes remove
  total_number_of_genes<-nrow(pseudo.object@assays$RNA$counts)
  numer_of_genes_rm<-total_number_of_genes-length(genes_to_keep)
  message(paste0("number of genes removed below cutoff: ",numer_of_genes_rm))
  return(subset(pseudo.object,features = genes_to_keep))
}

#---------function to normalize counts with DESeq2 function based on a seurat object of pseudocounts/metacells------------
NormalizeData_DESeq2<-function(sc.metacell,assay="DESeq2"){
  #get count data
  average.exp <- sc.metacell@assays[[assay]]$counts
  #create DESeq2 object
  dds <- DESeqDataSetFromMatrix(average.exp,
                                colData = data.frame(colnames(average.exp)),
                                design = ~1)
  #normalize with DESeq2
  dds<-DESeq(dds)
  #get norm counts and add to data slot of seurat metacell object
  normalized_counts <- counts(dds, normalized = TRUE)
  sc.metacell@assays[[assay]]@data<-normalized_counts
  
  return(sc.metacell)
}

#---------function to add vst transformed counts to the scale.data slot with DESeq2 function based on a seurat object of pseudocounts/metacells------------
#remove batch effects with limma http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
ScaleData_DESeq2<-function(sc.metacell,assay="DESeq2",conserve_covariate=NULL,batch1=NULL,batch2=NULL,regress_covariates=NULL){
  #comprehensive function to check if values in a vector are categorical or continues
  is_categorical <- function(vec) {
    # Check if the vector is a factor or character
    if (is.factor(vec) || is.character(vec)) {
      return(TRUE)
    }
    # If numeric, check the number of unique values
    if (is.numeric(vec)) {
      unique_values <- length(unique(vec))
      total_values <- length(vec)
      # Set a threshold for determining if it's categorical
      threshold <- 0.1 * total_values
      return(unique_values < threshold)
    }
    # Default to FALSE for other types
    return(FALSE)
  }
  
  message(paste0("returning vst transformed counts in scale data slot"))
  #get count data
  average.exp <- sc.metacell@assays[[assay]]@counts
  #create DESeq2 object
  if(!is.null(conserve_covariate)){
    if (length(conserve_covariate)!=1) {stop("conserve_covariate needs to be a single feature and should be the feature of interesst, like treatment or disease")    }
    message(paste0(conserve_covariate, "used as conserve_covariate in dds"))
    meta.data<-sc.metacell@meta.data
    dds <- DESeqDataSetFromMatrix(average.exp,
                                  colData = meta.data,
                                  design = as.formula(paste("~",conserve_covariate)))
  }else{
    message("no conserve_covariate used in dds, just intercept")
    dds <- DESeqDataSetFromMatrix(average.exp,
                                  colData = data.frame(colnames(average.exp)),
                                  design = ~1)
  }
  #perform variance stabilization
  average.exp.vst <- tryCatch({
    vst(dds, blind = TRUE)
  }, error = function(e) {
    message("vst() failed: ", e$message, "\nFalling back to varianceStabilizingTransformation()")
    varianceStabilizingTransformation(dds, blind = TRUE)
  })
  
  #in addition of correction by vst normalization limmaBatchCorretion can be used to remove categorical batch effects (1/2) and regress confounding covariates
  if (!is.null(batch1)|!is.null(regress_covariates)) {
    #get vst data and meta data
    mat <- assay(average.exp.vst)
    meta.data<-sc.metacell@meta.data
    
    #check if one of the variables to include in the model is not in meta.data, in that case check if its in the reductions
    if (!all(c(batch1,batch2,regress_covariates) %in% colnames(meta.data))) {
      message(paste0( c(batch1,batch2,regress_covariates)[!c(batch1,batch2,regress_covariates) %in% colnames(meta.data)]," not in meta.data \n"))
      for (reduction in names(sc.metacell@reductions)) {
        if (reduction==names(sc.metacell@reductions)[1]) {
          reductions_df<-Embeddings(sc.metacell,reduction = reduction)
        }else{
          reductions_df<-cbind(reductions_df,Embeddings(sc.metacell,reduction = reduction))
        }
      }
      #if the missing meta.data is in reductions, add all reduction to meta data  
      if (all(c(batch1,batch2,regress_covariates) %in% c(colnames(reductions_df),colnames(meta.data)))) {
        message(paste0( c(batch1,batch2,regress_covariates)[!c(batch1,batch2,regress_covariates) %in% colnames(meta.data)]," found in reductions \n"))
        meta.data<-cbind(meta.data,reductions_df)
      }else{
        message(paste0(c(batch1,batch2,regress_covariates)[!c(batch1,batch2,regress_covariates) %in% c(colnames(reductions_df),colnames(meta.data))]," also not in reductions, cant be found")) 
        stop()
      }
    }#end of if variables in meta.data
    
    #create design matrix indicating the relevant covariates to be preserved (conserve_covariate)
    message(paste0("preserving covariates: ",conserve_covariate))
    if (is.null(conserve_covariate)) {
      #only intercept
      message("no conserve_covariate, therefore model matrix with only the intercept")
      design <- model.matrix(~1, data = meta.data) 
    }else{
      design <- model.matrix(as.formula(paste("~1 + ",conserve_covariate)),meta.data)
    }
    
    #check if batch (1/2) needs to be removed and if its categorical
    if (!is.null(batch1)){
      if (is_categorical(batch1)) {
        message(paste0("running limma batch1 removal for: ",batch1))
        batch1<-meta.data[[batch1]]
      }else{
        message(paste0(batch1," is not categorical!!! cant be used as batch1"))
      }
    }
    if (!is.null(batch2)){
      if (is_categorical(batch2)) {
        message(paste0("running limma batch2 removal for: ",batch2))
        batch2<-meta.data[[batch2]]
      }else{
        message(paste0(batch2," is not categorical!!! cant be used as batch2"))
      }
    }
    
    #check if all features in regress_covariates are numeric values in a matrix, otherwise transform categorical to numeric
    if (!is.null(regress_covariates)) {
      message(paste0(c("running regression for: ",paste0(regress_covariates,", "))))
      for (feature in regress_covariates) {
        if (is_categorical(meta.data[[feature]])) {
          message(paste0(feature," looks categorical and is converted to numeric"))
          meta.data[feature][is.na(meta.data[feature])]<-"NA"
          meta.data[feature]<-as.numeric(as.factor(meta.data[[feature]]))
        }
      }
      regress_covariates<-meta.data[regress_covariates]
    }
    
    #run correction
    mat <- removeBatchEffect(mat, batch = batch1 , batch2 = batch2 , covariates = regress_covariates , design = design)
    assay(average.exp.vst) <- mat
    
    #end of if limmaBatchCorrection is run 
  }else{
    message("no regression")
  }
  
  #add vst to scale data slot
  sc.metacell@assays[[assay]]@scale.data<-assay(average.exp.vst)
  return(sc.metacell)
}

#---------function to find highly variable features with DESeq2 function based on a seurat object of pseudocounts/metacells------------
FindVariableFeatures_DESeq2<-function(sc.metacell,nHVG=500,assay="DESeq2"){
  message(paste0("calculating the ",nHVG," highly variable genes after vst"))
  #get scaled count data
  average.exp.vst <- sc.metacell@assays[[assay]]@scale.data
  # #create DESeq2 object
  # dds <- DESeqDataSetFromMatrix(average.exp,
  #                               colData = data.frame(colnames(average.exp)),
  #                               design = ~1)
  # #perform variance stabilization
  # average.exp.vst <- vst(dds, blind=TRUE)
  #get n highly variable genes
  rv <- matrixStats::rowVars(average.exp.vst,useNames = F)
  select <- order(rv, decreasing=TRUE)[seq_len(min(nHVG, length(rv)))]
  HVG <- rownames(average.exp.vst[select,])
  #set HVG in the object
  sc.metacell@assays[[assay]]@var.features<-HVG
  return(sc.metacell)
}

#---------function to calculate PCA based on vst with DESeq2 and return as reduction in seurat object of pseudocounts/metacells------------
RunPCA_DESeq2<-function(sc.metacell,assay="DESeq2"){
  message(paste0("calculate PCA based on variable features slot and scale.data"))
  #get variable features
  var_feat<-VariableFeatures(sc.metacell,assay = assay)
  scaled_count.data<-sc.metacell@assays[[assay]]@scale.data[var_feat,]
  #calculate PCs
  message(paste0("number of var feat used: ",length(var_feat)))
  pca <- prcomp(t(sc.metacell@assays[[assay]]@scale.data[var_feat,]))
  #add DESeq2 pca to reduction
  sc.metacell[[paste0(assay,"_pca")]] <- CreateDimReducObject(embeddings = pca$x, loadings = pca$rotation,stdev = pca$sdev, key = paste0(assay,"_PCA_"), assay = assay)
  return(sc.metacell)
}


#-------ratio of celltype composition plot with replicates-------------
composition.ratio.box.plot <- function(sc.object,group1="seurat_clusters",group2,sampleIDs="orig.ident",cols=safe_c2,exclude_ct=NULL,capY=NULL,label_outlier=F,outlier_cutoff=3,
                                       add_statistic=T,print_stats = F,normalize.by="total",p_val_style="p.adj.signif",p.cut_off=0.05,return_df=F){
  #create a meta data slot were the condition label (usually group2) gets split into the underlying replicates (samples)
  sc.object$tmp_group2_sampleID<-paste0(unlist(sc.object[[group2]]),"_",unlist(sc.object[[sampleIDs]]))
  #get the cell numbers
  breakdown<-table(sc.object[[c(group1,"tmp_group2_sampleID")]])
  #normalize each sample to equal input
  #normalize.by defines what to use as reference. usually the total cell number of a sample (default)
  #can be changed to a reference
  if (normalize.by=="total") {
    breakdown<-apply(breakdown, 2, function(x){(x/sum(x))*100 })
  }else{
    #normalize the cell state composition by the total number of the higher level celltype annotation
    if (normalize.by=="celltypes.short") {
      #get ref list of celltype and cellstate relation
      df_cell_anno <- sc.object@meta.data[c(group1,"celltypes.short")]
      df_cell_anno <- split(as.character(df_cell_anno[[group1]]), df_cell_anno$celltypes.short)
      df_cell_anno <- lapply(df_cell_anno, unique)
      #correct EndoEC and LymEC
      df_cell_anno$EndoEC<-c("EndoEC")
      df_cell_anno$LYMEC<-c("LymEC")
      df_cell_anno$EC<-c("capillaryEC", "angioEC",  "veinEC" , "arteryEC" )
      #isolate relevant number and calculate freq
      for (ct in names(df_cell_anno)) {
        if (ct==names(df_cell_anno)[1]) {
          df <- breakdown[rownames(breakdown) %in% df_cell_anno[[ct]],]
          combined_df<-apply(df, 2, function(x){(x/sum(x))*100 })
        }else{
          df <- breakdown[rownames(breakdown) %in% df_cell_anno[[ct]],]
          if (is.null(nrow(df))) {
            message(paste0("cell state: ",ct," is the only state of that type, therefore normalized to total cell number"))
            df<-t(as.data.frame((df / colSums(breakdown))*100))
            rownames(df) <- df_cell_anno[[ct]]
          } else {
            df<-apply(df, 2, function(x){(x/sum(x))*100 })
          }
          combined_df<-rbind(combined_df,df)
        }
      }
      breakdown<-combined_df
    }else{stop("normalize.by not impremented for this")}
  }
  breakdown.df = as.data.frame.matrix(breakdown)
  breakdown.df = melt(t(breakdown.df))
  names(breakdown.df)<-c("tmp_group2_sampleID","group1","cell_number_pt")
  breakdown.df$group2<-as.character(breakdown.df$tmp_group2_sampleID)
  breakdown.df$group1<-as.character(breakdown.df$group1)
  #add column with only condition lable
  for (type_group2 in unique(unlist(sc.object[[group2]]))) {
    breakdown.df$group2[str_detect(breakdown.df$group2,type_group2)]<-type_group2
  }
  #add column with only the patient ID (to be able to find pairs)
  breakdown.df$sample<-sapply( str_split(breakdown.df$tmp_group2_sampleID,"_"), function(x) return(x[3]))
  #check which patients have pairs in the sc data
  samples<-as.character(unlist(unique(sc.object[[sampleIDs]])))
  samples<-gsub("_PRE","",samples)
  samples<-gsub("_POST","",samples)
  samples_with_pairs<-names(table(samples))[table(samples)>1]
  #reduce cell number table to patients that still have pairs
  breakdown.df<-breakdown.df[breakdown.df$sample %in% samples_with_pairs,]
  #add column with only the time point (pre or post )
  breakdown.df$Tp<-sapply(str_split(breakdown.df$tmp_group2_sampleID,"_"), function(x) return(x[4]))
  #add colmun with only the study group (IVIG or PLACEO)
  breakdown.df$treament<-breakdown.df$group2
  breakdown.df$treament<-gsub("_POST","", breakdown.df$treament)
  breakdown.df$treament<-gsub("_PRE","", breakdown.df$treament)
  if(return_df){return(breakdown.df)}
  #calculate the ratio for each patient, per cell type
  ratios<-NULL
  for (celltype in unique(breakdown.df$group1)) {
    for (ID in unique(breakdown.df$sample)) {
      T0<-breakdown.df$cell_number_pt[breakdown.df$group1==celltype & breakdown.df$sample==ID &breakdown.df$Tp=="0B"]
      T1<-breakdown.df$cell_number_pt[breakdown.df$group1==celltype & breakdown.df$sample==ID &breakdown.df$Tp=="1B"]
      ratios[[celltype]][[ID]]<-T1/T0
    }
  }
  ratios<-melt(ratios)
  #remove Inf values where the celltype did not exist in one of the biopsies 
  ratios<-ratios[ratios$value!="Inf",]
  
  #!! IMPORTANT transform the ratios to log(ratio) to scale both directions of change equaly!!
  ratios$value<-log(ratios$value)
  #rm now where log(0) is -Inf
  ratios<-ratios[ratios$value!="-Inf",]
  
  #add the study group back based on the sampleID
  ratios$Condition<-ratios$L2
  for (ID in unique(ratios$L2)) {
    ratios$Condition[ratios$L2==ID]<-breakdown.df$treament[breakdown.df$sample==ID]
  }
  #give the ct (L1) their order as factor
  if (is.factor(sc.object[[c("tmp_group2_sampleID",group1)]][[2]])) {
    ratios$L1<-factor(ratios$L1,levels = levels(sc.object[[c("tmp_group2_sampleID",group1)]][,2]))
  }else{
    ratios$L1<-factor(ratios$L1)
  }
  
  #exclue celltypes as required (optional input)
  if (!is.null(exclude_ct)) {
    ratios<-ratios[!ratios$L1 %in% exclude_ct,]
    ratios$L1<-factor(ratios$L1)}
  #add order to study group
  ratios$Condition<-factor(ratios$Condition,levels = c("PLACEBO","IVIG"))
  
  #Calculate extrem outliers (Q1, Q3, IQR *3) outlier_cutoff
  ratios_with_outliers <- ratios[!is.na(ratios$value),] %>%
    group_by(L1) %>%  # Assuming L1 is the grouping factor
    mutate(Q1 = quantile(value, 0.25),
           Q3 = quantile(value, 0.75),
           IQR = Q3 - Q1,
           lower_bound = Q1 - outlier_cutoff * IQR,
           upper_bound = Q3 + outlier_cutoff * IQR,
           is_outlier = ifelse(value < lower_bound | value > upper_bound, TRUE, FALSE)) %>% # Flag outliers
    as.data.frame()
  #add column to label outlier
  ratios_with_outliers$is_outlier_label<-""
  ratios_with_outliers$is_outlier_label[ratios_with_outliers$is_outlier]<-ratios_with_outliers$L2[ratios_with_outliers$is_outlier]
  
  #(optional input) sometime there are extrem ratio values, in that case exclude them by setting y max to 10
  if (!is.null(capY)) { 
    y_lim_manual<-c(-capY, capY)
    y_break_manual <-seq(-capY, capY, by = 1)
    message(paste0("number of ratios > 10: ",sum(na.omit(ratios$value) > capY)))
  }else{
    y_lim_manual<-c(round(min(na.omit(ratios$value)),0), round(max(na.omit(ratios$value)),0))
    y_break_manual <-seq(round(min(na.omit(ratios$value)),0), round(max(na.omit(ratios$value)),0), by = 1)
  }
  
  #label outliers by sample ID (optional)
  if (label_outlier) {
    p<- ggplot(ratios_with_outliers, aes(x=L1, y=value)) + 
      geom_boxplot(aes(fill=Condition),position=position_dodge2(padding=0.1), pch=21, show.legend = T,outlier.shape = NA)+
      geom_point(aes(fill=Condition),position=position_jitterdodge(jitter.width=0, dodge.width = 0.8),pch=21,  show.legend = T)+
      scale_fill_manual(values =cols)+
      geom_text_repel(aes(label=is_outlier_label,fill=Condition),
                      position=position_jitterdodge(jitter.width=0, dodge.width = 0.8),
                      size=3,min.segment.length = 0) +  
      ylab("logFC of pct cells per sample")+
      geom_hline(yintercept = 0,linetype="dashed")+
      scale_y_continuous(limits = y_lim_manual, breaks = y_break_manual)+
      theme_classic()+
      theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45,hjust = 1),
            axis.text = element_text(size=15),
            axis.title = element_text(size=15))
  }else{
    p<-  ggplot(ratios_with_outliers, aes(x=L1, y=value)) + 
      geom_boxplot(aes(fill=Condition),position=position_dodge2(padding=0.1), pch=21, show.legend = T,outlier.shape = NA)+
      geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.8),pch=21, aes(fill=Condition), show.legend = T)+
      scale_fill_manual(values =cols)+
      ylab("logFC of pct cells per sample")+
      geom_hline(yintercept = 0,linetype="dashed")+
      scale_y_continuous(limits = y_lim_manual, breaks = y_break_manual)+
      theme_classic()+
      theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45,hjust = 1),
            axis.text = element_text(size=15),
            axis.title = element_text(size=15))
  }
  #if normalization by CM, overwrite yaxis lable
  if (normalize.by=="celltypes.short") { p<- p+ylab("cell type logFC normalized to celltype per sample")}
  
  
  #add statistics (optional)
  if(add_statistic){
    library(rstatix)
    # Statistical test
    # Perform pairwise comparisons using Wilcoxon tests
    stat_test <- ratios_with_outliers[!ratios_with_outliers$is_outlier,]  %>%
      group_by(L1) %>%
      pairwise_wilcox_test(value ~ Condition, p.adjust.method = "BH")  %>% 
      add_significance("p.adj")  %>%
      add_xy_position(x = "L1",step.increase = 0.05)
    #add column with only significant pvalues
    stat_test$p.adj.signif.num<-stat_test$p.adj
    stat_test$p.adj.signif.num[stat_test$p.adj.signif.num > p.cut_off]<-"ns"
    stat_test$p.adj.signif.NA<-stat_test$p.adj
    stat_test$p.adj.signif.NA[stat_test$p.adj.signif.NA > p.cut_off]<-NA
    #show test results (optional)
    if (print_stats) {View(stat_test)}
    #add to plot
    p <- p + stat_pvalue_manual(stat_test, label = p_val_style, tip.length = 0.01) #p.adj.signif, p.adj or p.adj.signif.num, p.adj.signif.NA
  }
  print(p)
}


#function to perform DGE to find marker genes of clusters in a metacell seurat object (based on limma voom)
FindAllMarkers_pseudobulk <-function(metacell.object,cluster_identifier="seurat_clusters",filter_pct_exp_per_cluster=F,filter_pct_exp_per_cluster_cutoff=0.3){
  require(limma)
  
  #check if order(levels)exist
  if(is.factor(FetchData(metacell.object,vars = cluster_identifier)[,1])){
    clusters<-levels(metacell.object[[cluster_identifier]][1,])
  }else{
    clusters<-unique(unlist(metacell.object[[cluster_identifier]]))
  }
  
  marker_genes_result<-NULL
  #for each cluster group compare against all other groups to find cluster specific marker genes
  for (cluster in clusters) {
    #duplicate
    metacell.object2<-metacell.object
    #get count matrix of pseudobulked samples 
    count.matrix<-GetAssayData(metacell.object2,slot = "counts",assay = "RNA")
    
    #optional filter for genes that are robustly expressed in the target cluster
    if (filter_pct_exp_per_cluster) {
      #subset to cluster
      Idents(metacell.object2)<-cluster_identifier
      metacell.object2_cluster<-subset(metacell.object2,ident=cluster)
      #only continue if there is sufficient samples pseudobulked
      if (ncol(metacell.object2_cluster)>1) {
        #check how often a gene has non-zero expression  
        how_often_in_sample <- rowSums(!metacell.object2_cluster@assays$RNA$counts == 0)
        names(how_often_in_sample) <- rownames(metacell.object2_cluster@assays$RNA$counts)
        how_often_in_sample_pct<-how_often_in_sample/ncol(metacell.object2_cluster@assays$RNA$counts)
        #genes below cutoff
        below_cutoff<-names(how_often_in_sample_pct)[how_often_in_sample_pct<filter_pct_exp_per_cluster_cutoff]
        message(paste0("for cluster: ",cluster,". ",length(below_cutoff)," genes are below cutoff ",filter_pct_exp_per_cluster_cutoff," for coverage and are rm for marker DE analysis"))
        #remove from count.matrix
        count.matrix<-count.matrix[!rownames(count.matrix)%in%below_cutoff,]
      }else{message(paste0("not enough cells skipped cluster :",cluster))}
    }
    
    
    #generate group vector 
    group<-FetchData(metacell.object2,vars = cluster_identifier)[,1]
    group<-as.character(group)
    #classify each samples as cluster of interesst or pool of all others
    group[group!=cluster]<-"other_clusters"
    group[group==cluster]<-"cluster_of_interesst"
    
    #variance control with voom 
    mm <- model.matrix(~0 + group)
    y <- voom(count.matrix, mm, plot = F)
    fit <- lmFit(y, mm)
    #DEGs
    contr <- makeContrasts(groupcluster_of_interesst - groupother_clusters, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    DEGs <- topTable(tmp, number = length(count.matrix[,1]), adjust.method = "fdr")
    #add genes as separat column
    DEGs$gene<-rownames(DEGs)
    #add which cluster the results are for
    DEGs$cluster<-cluster
    #save overall result
    marker_genes_result<-rbind(marker_genes_result,DEGs)
  }
  return(marker_genes_result)
}



