#' CMSPlus main function
#'
#' @param exp2symbol exp2symbol a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols.
#' @param prob A numeric value between 0 and 1 specifying the minimum posterior probability threshold. Required to assign a predicted subtype. Samples with maximum subtype probability below this threshold will be classified as "unassigned". Default is 0.5.
#' @param plot TRUE produces plots; FALSE suppresses plotting.
#' @param CMSlabel A character string specifying the label used in plots to denote CMS subtypes.
#' @param CMSPluslabel A character string specifying the label used in plots to denote CMSPlus subtypes.
#' @param parallel.sz Integer specifying the number of parallel workers used for GSVA computation. Default is 1.
#' @param InterGroupRandomize Logical value indicating whether to perform randomization in plots.
#' @param seed Integer specifying the random seed used to ensure reproducibility of randomization step.
#'
#' @return list including CMSPlusLabel, CMSLabel, gsva_matrix and HeatmapPlot
#' @export
#'
#' @import GSVA
#' @import CMSclassifier
#' @import ComplexHeatmap
#' @import circlize
#' @import parallel
#' @import randomForest
#' @import ranger
CMSPlus <- function(exp2symbol, prob=0.5, plot=TRUE, CMSlabel="RF.nearestCMS", CMSPluslabel="nearest", parallel.sz=0, InterGroupRandomize=T, seed=123){
  message('Performing GSVA......')
  gsva_matrix <- gsva(as.matrix(exp2symbol), CMSPlusSignatures, method=c("gsva"), kcdf = "Gaussian", min.sz = 1, parallel.sz=parallel.sz)
  message('Done GSVA')
  gsva_matrix_t <- as.data.frame(t(gsva_matrix))
  colnames(gsva_matrix_t) <- gsub("-", "_", colnames(gsva_matrix_t))
  colnames(gsva_matrix_t) <- gsub(" ", "_", colnames(gsva_matrix_t))
  message('CMSPlus prediction......')
  cmsplus.prob = as.data.frame(predict(CMSPlusModel, data = gsva_matrix_t)$predictions)
  rownames(cmsplus.prob) <- colnames(exp2symbol)
  cmsplus.prob$nearest <- apply(cmsplus.prob[, 1:5], 1, function(x) {
    labs <- names(x)[x == max(x)]
    paste(labs, collapse = ",")
  })
  cmsplus.prob$predict <- apply(cmsplus.prob[, 1:5], 1, function(x) {
    labs <- names(x)[x > prob]
    if (length(labs) == 0) {
      "mix"
    } else {
      paste(labs, collapse = ",")
    }
  })
  message('CMSPlus prediction done')
  exp2symbol_entrez <- ChangeEntrezFromName(GeneCount=exp2symbol)
  Rfcms <- CMSclassifier::classifyCMS(exp2symbol_entrez,method="RF")[[3]]
  phenotype_df <- cbind(cmsplus.prob, Rfcms)
  if(plot){
    phenotype_df2 <- phenotype_df
    phenotype_df2[, CMSPluslabel] <- sapply(strsplit(phenotype_df2[, CMSPluslabel], "\\,"), function(x) x[1])
    phenotype_df2[, CMSlabel] <- sapply(strsplit(phenotype_df2[, CMSlabel], "\\,"), function(x) x[1])
    HeatmapPlot <- PlotHeatMap(ClinData=phenotype_df2, GSVAData=gsva_matrix, CMSPlus=CMSPluslabel, CMS=CMSlabel, SDI=NULL, InterGroupRandomize=InterGroupRandomize, seed=seed)
    message('Done Heatmap')
    HeatmapPlot
    res <- list(CMSPlusLabel=cmsplus.prob, CMSLabel=Rfcms, gsva_matrix=gsva_matrix, HeatmapPlot=HeatmapPlot)
  } else {
    res <- list(CMSPlusLabel=cmsplus.prob, CMSLabel=Rfcms, gsva_matrix=gsva_matrix)
  }
  return(res)
}
#'  Changing the expression profile from Gene symbol to Entrez ID
#'
#' @param GeneCount a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols.
#'
#' @return a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to Entrez ID
#' @export
#'
ChangeEntrezFromName <- function(GeneCount){
  GeneNameb1 <- GeneCount[as.character(intersect(rownames(GeneCount), Ensembl2Name2NCBI$Gene.name[which(!is.na(Ensembl2Name2NCBI$NCBI.gene.ID))])), ]
  GeneNameb1Rowname <- as.character(Ensembl2Name2NCBI$NCBI.gene.ID[match(rownames(GeneNameb1), Ensembl2Name2NCBI$Gene.name)])
  Convert <- DeleteDupChoseMax(GeneNameb1, GeneNameb1Rowname)
  return(Convert)
}
#' Deleting duplication genes
#'
#' @param ExpDataframe a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols.
#' @param ReplaceName a character with entrez id to be replaced.
#' @param by delete duplications by column or row.
#' @param ori keep original gene symbol.
#'
#' @return a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to entrez ID.
#' @export
#'
DeleteDupChoseMax <- function(ExpDataframe, ReplaceName, by="row", ori=FALSE){
  require(magrittr)
  if(by == "row"){
    ReplaceRowname <- ReplaceName
    dupchar <- ReplaceRowname[duplicated(ReplaceRowname)] %>% unique
    if(length(dupchar) > 0){
      expUnique <- ExpDataframe[!ReplaceRowname %in% dupchar, ]
      rownames(expUnique) <- ReplaceRowname[!ReplaceRowname %in% dupchar]
      expDup <- do.call(rbind, lapply(dupchar, function(x){
        Single <- ExpDataframe[which(ReplaceRowname == x), ]
        tmp <- Single[which.max(rowMeans(Single)), ]
        return(tmp)
      }))
      rownames(expDup) <- dupchar
      expDup <- RowDFtoRawNum(expDup)
      FinalDataframe <- rbind(expUnique, expDup)
    } else {
      FinalDataframe <- ExpDataframe
      rownames(FinalDataframe) <- ReplaceRowname
    }
  }
  if(by == "col"){
    ReplaceColname <- ReplaceName
    dupchar <- ReplaceColname[duplicated(ReplaceColname)] %>% unique
    expUnique <- ExpDataframe[, !ReplaceColname %in% dupchar]
    colnames(expUnique) <- ReplaceColname[!ReplaceColname %in% dupchar]
    expDup <- do.call(cbind, lapply(dupchar, function(x){
      Single <- ExpDataframe[, which(ReplaceColname == x)]
      tmp <- Single[, which.max(colMeans(Single))]
      return(tmp)
    }))
    colnames(expDup) <- dupchar
    FinalDataframe <- cbind(expUnique, expDup)
    if(ori == TRUE){
      uniquenames <- colnames(ExpDataframe)[!ReplaceColname %in% dupchar]
      remains <- do.call(c, lapply(dupchar, function(x){
        Single <- ExpDataframe[, which(ReplaceColname == x)]
        tmp <- colnames(Single)[which.max(colMeans(Single))]
        return(tmp)
      }))
      colnames(FinalDataframe) <- c(uniquenames, remains)
    }
  }
  return(FinalDataframe)
}
#' row to row
#'
#' @param ExpDataframe a dataframe with Gene Expression Profiles data values.
#'
#' @return a dataframe with Gene Expression Profiles data values
#' @export
#'
RowDFtoRawNum <- function(ExpDataframe){
  for(i in 1:dim(ExpDataframe)[2]){
    tmp <- as.numeric(ExpDataframe[, i])
    if(i == 1){
      Matrix <- tmp
    } else {
      Matrix <- cbind(Matrix, tmp)
    }
  }
  colnames(Matrix) <- colnames(ExpDataframe)
  rownames(Matrix) <- rownames(ExpDataframe)
  df <- as.data.frame(Matrix)
  return(df)
}
#' Ploting heatmap
#'
#' @param ClinData a data frame with CMSPlus labels and CMS labels.
#' @param GSVAData GSVA matrix
#' @param CMSPlus CMSPlus label name
#' @param CMS CMS label name
#' @param SDI SDI label name
#' @param InterGroupRandomize TRUE or FALSE
#' @param seed numeric
#'
#' @return list including CMSPlusLabel, CMSLabel, gsva_matrix and HeatmapPlot
#' @export
#'
#' @import ComplexHeatmap
#' @import circlize
PlotHeatMap <- function(ClinData, GSVAData, CMSPlus, CMS, SDI, InterGroupRandomize=T, seed=123){
  columnAnno = HeatmapAnnotation(
    CMSPlus = ClinData[, CMSPlus],
    CMS=ClinData[, CMS],
    col = list(CMSPlus =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-TME+"='#3C5488', "CMS4-TME-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"))
  )
  col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
  col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))
  dendtcga = cluster_within_group(GSVAData, ClinData[, CMSPlus])
  P1 <- Heatmap(GSVAData, col = col_fun, cluster_rows = FALSE, cluster_columns = dendtcga, top_annotation = columnAnno, column_labels=rep("", dim(GSVAData)[2]))
  if(InterGroupRandomize){
    set.seed(seed)
    columnOrder <- c(sample(intersect(column_order(P1), grep("CMS1", ClinData[, CMSPlus]))), sample(intersect(column_order(P1), grep("CMS2", ClinData[, CMSPlus]))), sample(intersect(column_order(P1), grep("CMS3", ClinData[, CMSPlus]))), sample(intersect(column_order(P1), grep("CMS4-TME\\+", ClinData[, CMSPlus]))), sample(intersect(column_order(P1), grep("CMS4-TME-", ClinData[, CMSPlus]))))
  } else {
    columnOrder <- c(intersect(column_order(P1), grep("CMS1", ClinData[, CMSPlus])), intersect(column_order(P1), grep("CMS2", ClinData[, CMSPlus])), intersect(column_order(P1), grep("CMS3", ClinData[, CMSPlus])), intersect(column_order(P1), grep("CMS4-TME\\+", ClinData[, CMSPlus])), intersect(column_order(P1), grep("CMS4-TME-", ClinData[, CMSPlus])))
  }
  if(length(SDI)>0){
    ClinData[, SDI] <- gsub("-", "_", ClinData[, SDI])
    columnAnno2 = HeatmapAnnotation(
      CMSPlus = ClinData[, CMSPlus][columnOrder],
      CMS=ClinData[, CMS][columnOrder],
      CMSSDI=ClinData[, SDI][columnOrder],
      col = list(CMSPlus =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-TME+"='#3C5488', "CMS4-TME-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"), CMSSDI=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighSDI='#C2CB1E', CMS4_LowSDI='#537A34')
      ))
  } else {
    columnAnno2 = HeatmapAnnotation(
      CMSPlus = ClinData[, CMSPlus][columnOrder],
      CMS=ClinData[, CMS][columnOrder],
      #CMSSDI=ClinData[, SDI][columnOrder],
      col = list(CMSPlus =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-TME+"='#3C5488', "CMS4-TME-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"), CMSSDI=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighSDI='#C2CB1E', CMS4_LowSDI='#537A34')
      ))
  }
  GSVAData2 <- GSVAData[c(1:4, 6:7, 5, 8:15, match(TMIGOrder[c(1:6, 9:27, 7:8), 1], rownames(GSVAData))), columnOrder]
  row.subsections1 <- c(3, 4, 5, 3)
  row_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4"), row.subsections1))
  col.subsections1 <- as.numeric(table(ClinData[, CMSPlus])[c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-")])
  col_split1 = data.frame(rep(names(table(ClinData[, CMSPlus])), col.subsections1))
  P2 <- Heatmap(GSVAData2[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = columnAnno2, column_labels=rep("", dim(GSVAData2)[2]), row_split = row_split1, column_split = col_split1)
  row.subsections2 <- c(6, 21)
  row_split2 = data.frame(factor(rep(c("Fibroblastes", "Immune"), row.subsections2), levels=c("Immune", "Fibroblastes")))
  P3 <- Heatmap(GSVAData2[16:42, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, column_labels=rep("", dim(GSVAData2)[2]), row_split = row_split2, column_split = col_split1)
  ht_list = P2%v%P3
  draw(ht_list)
}

