## 添加了绘制点颜色的功能，相比于原来的特点。
tcgaCompare_test <- function (maf, capture_size = NULL, tcga_capture_size = 35.8, 
          cohortName = NULL, tcga_cohorts = NULL, primarySite = FALSE, 
          col = c("gray70", "black"), bg_col = c("#EDF8B1", "#2C7FB8"), 
          medianCol = "red", decreasing = FALSE, logscale = TRUE, rm_hyper = FALSE, 
          rm_zero = TRUE, cohortFontSize = 0.8, axisFontSize = 0.8) 
{
  # test.R 文件生成了tcga_cohort.txt.gz这个文件
  # 这个文件没有细节
  # 缩略图如下
  # cohort  Tumor_Sample_Barcode    total   site
  # ACC     TCGA-OR-A5KB-01A-11D-A30A-10    1803    Adrenal Gland Carcinoma
  # ACC     TCGA-PK-A5HB-01A-11D-A29I-10    962     Adrenal Gland Carcinoma
  # ACC     TCGA-OR-A5JA-01A-11D-A29I-10    464     Adrenal Gland Carcinoma
  # ACC     TCGA-OR-A5LJ-01A-11D-A29I-10    395     Adrenal Gland Carcinoma
  
  isac_meta <- data.frame(getSampleSummary(maf_merge)$Tumor_Sample_Barcode, rep(1, length(getSampleSummary(maf_merge)$Tumor_Sample_Barcode)))
  colnames(isac_meta) <- c("Tumor_Sample_Barcode", "age_year")
  tcga.cohort = system.file("extdata", "tcga_cohort.txt.gz", 
                            package = "maftools")
  tcga.cohort = data.table::fread(file = tcga.cohort, sep = "\t", 
                                  stringsAsFactors = FALSE)
  table_age <- read.table("~/OneDrive - Karolinska Institutet/personal_files/PhD/Project/ISAC/results/table/age_tcga_cancer_2023-09-05_21-59-01.067349.csv", header = T)
  table_age <- table_age[,-1]
  tcga.cohort <- merge(tcga.cohort, table_age, by = "Tumor_Sample_Barcode")
  if (primarySite) {
    tcga.cohort = tcga.cohort[, .(Tumor_Sample_Barcode, total, 
                                  site)]
    colnames(tcga.cohort)[3] = "cohort"
  }
  else {
    tcga.cohort = tcga.cohort[, .(Tumor_Sample_Barcode, total, 
                                  cohort, age)]
  }
  if (!is.null(tcga_cohorts)) {
    tcga.cohort = tcga.cohort[cohort %in% tcga_cohorts]
    if (nrow(tcga.cohort) == 0) {
      stop("Something went wrong. Provide correct names for 'tcga_cohorts' arguments")
    }
  }
  if (length(maf) == 1) {
    maf = list(maf)
  }
  maf.mutload = lapply(maf, function(m) {
    x = getSampleSummary(m)[, .(Tumor_Sample_Barcode, total)]
    if (rm_zero) {
      warning(paste0("Removed ", nrow(x[x$total == 0]), 
                     " samples with zero mutations."))
      x = x[!total == 0]
    }
    x
  })
  if (is.null(cohortName)) {
    cohortName = paste0("Input", seq_len(length(maf)))
  }
  else if (length(cohortName) != length(maf)) {
    stop("Please provide names for all input cohorts")
  }
  names(maf.mutload) = cohortName
  # browser()
  maf.mutload = data.table::rbindlist(l = maf.mutload, idcol = "cohort")
  tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
  maf.mutload$total = as.numeric(as.character(maf.mutload$total))
  if (rm_hyper) {
    tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
    tcga.cohort = lapply(tcga.cohort, function(x) {
      xout = boxplot.stats(x = x$total)$out
      if (length(xout) > 0) {
        message(paste0("Removed ", length(xout), " outliers from ", 
                       x[1, cohort]))
        x = x[!total %in% xout]
      }
      x
    })
    tcga.cohort = data.table::rbindlist(l = tcga.cohort)
    xout = boxplot.stats(x = maf.mutload$total)$out
    if (length(xout) > 0) {
      message(paste0("Removed ", length(xout), " outliers from Input MAF"))
      maf.mutload = maf.mutload[!total %in% xout]
    }
  }
  if (is.null(capture_size)) {
    tcga.cohort = rbind(tcga.cohort, maf.mutload)
    pt.test = pairwise.t.test(x = tcga.cohort$total, g = tcga.cohort$cohort, 
                              p.adjust.method = "fdr")
    pt.test.pval = as.data.frame(pt.test$p.value)
    data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
    colnames(pt.test.pval)[1] = "Cohort"
    message("Performing pairwise t-test for differences in mutation burden..")
    pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
    colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
    tcga.cohort[, `:=`(plot_total, total)]
  }
  else {
    message(paste0("Capture size [TCGA]:  ", tcga_capture_size))
    message(paste0("Capture size [Input]: ", capture_size))
    maf.mutload[, `:=`(total_perMB, total/capture_size)]
    tcga.cohort[, `:=`(total_perMB, total/tcga_capture_size)]
    
    ## 这里需要给isac的项目也添加上年龄信息
    library(tidyverse)
    # isac_meta <- read.table("~/OneDrive - Karolinska Institutet/personal_files/PhD/Project/ISAC/results/table/tmb_merge_with_meta_data_shared_samples_2023-08-09_14-19-53.08741.csv", 
    #                         header = T, sep = ",")     %>%
    #  dplyr::select(Tumor_Sample_Barcode, age_year)
    
    
    maf.mutload <- merge(maf.mutload, isac_meta, by = "Tumor_Sample_Barcode")
    colnames(maf.mutload)[5] <- "age"
    print("----------------------")
    tcga.cohort = rbind(tcga.cohort, maf.mutload)
    message("Performing pairwise t-test for differences in mutation burden (per MB)..")
    pt.test = pairwise.t.test(x = tcga.cohort$total_perMB, 
                              g = tcga.cohort$cohort, p.adjust.method = "fdr")
    pt.test.pval = as.data.frame(pt.test$p.value)
    data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
    colnames(pt.test.pval)[1] = "Cohort"
    pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
    colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
    tcga.cohort[, `:=`(plot_total, total_perMB)]
  }
  tcga.cohort.med = tcga.cohort[, .(.N, median(plot_total)), 
                                cohort][order(V2, decreasing = decreasing)]
  tcga.cohort$cohort = factor(x = tcga.cohort$cohort, levels = tcga.cohort.med$cohort)
  colnames(tcga.cohort.med) = c("Cohort", "Cohort_Size", "Median_Mutations")
  tcga.cohort$TCGA = ifelse(test = tcga.cohort$cohort %in% 
                              cohortName, yes = "Input", no = "TCGA")
  tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
  plot.dat = lapply(seq_len(length(tcga.cohort)), function(i) {
    x = tcga.cohort[[i]]
    x = data.table::data.table(rev(seq(i - 1, i, length.out = nrow(x))), 
                               x[order(plot_total, decreasing = TRUE), plot_total], 
                               x[, TCGA],
                               x[order(plot_total, decreasing = TRUE), age])
    x
  })
  names(plot.dat) = names(tcga.cohort)
  #browser()
  if (logscale) {
    y_lims = range(log10(data.table::rbindlist(l = plot.dat)[V2 != 
                                                               0][, V2]))
  }
  else {
    y_lims = range(data.table::rbindlist(l = plot.dat)[, 
                                                       V2])
  }
  y_min = floor(min(y_lims))
  y_max = ceiling(max(y_lims))
  y_lims = c(y_min, y_max)
  y_at = pretty(y_lims)
  # browser()
  plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, 
       axes = FALSE, xlab = NA, ylab = NA, frame.plot = TRUE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = grDevices::adjustcolor(col = "gray", alpha.f = 0.1))
  rect(xleft = seq(0, length(plot.dat) - 1, 1), ybottom = min(y_lims), 
       xright = seq(1, length(plot.dat), 1), ytop = y_max, col = grDevices::adjustcolor(col = bg_col, 
                                                                                        alpha.f = 0.2), border = NA)
  abline(h = pretty(y_lims), lty = 2, col = "gray70")
  #browser()
  lapply(seq_len(length(plot.dat)), function(i) {
    x = plot.dat[[i]]
    if (x[1, V3] == "Input") {
      if (logscale) {
        points(x$V1, log10(x$V2), pch = 16, cex = 0.4, 
               col = col[2])
      }
      else {
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[2])
      }
    }
    else {
      if (logscale) {
        # 先绘制非红色的点
        non_red <- x$V4 >= 18
        points(x$V1[non_red], log10(x$V2[non_red]), pch = 16, cex = 0.4, col = col[1])
        
        # 再绘制红色的点
        red <- x$V4 < 18
        points(x$V1[red], log10(x$V2[red]), pch = 16, cex = 0.4, col = "black")
      }
      else {
        points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[1])
      }
    }
  })
  samp_sizes = lapply(plot.dat, nrow)
  axis(side = 1, at = seq(0.5, length(plot.dat) - 0.5, 1), 
       labels = names(plot.dat), las = 2, tick = FALSE, line = -0.8, 
       cex.axis = cohortFontSize)
  axis(side = 3, at = seq(0.5, length(plot.dat) - 0.5, 1), 
       labels = unlist(samp_sizes), las = 2, tick = FALSE, line = -0.8, 
       font = 3, cex.axis = cohortFontSize)
  tcga.cohort.med[, `:=`(Median_Mutations_log10, log10(Median_Mutations))]
  if (decreasing) {
    sidePos = 2
    linePos = 2
  }
  else {
    sidePos = 2
    linePos = 2
  }
  if (logscale) {
    if (is.null(capture_size)) {
      axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, 
           labels = 10^(y_at), cex.axis = axisFontSize)
      mtext(text = "TMB", side = sidePos, line = linePos)
    }
    else {
      axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, 
           labels = 10^(y_at), cex.axis = axisFontSize)
      mtext(text = "TMB (per MB)", side = sidePos, line = linePos)
    }
  }
  else {
    if (is.null(capture_size)) {
      axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, 
           cex.axis = axisFontSize)
      mtext(text = "TMB", side = sidePos, line = linePos)
    }
    else {
      axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, 
           cex.axis = axisFontSize)
      mtext(text = "TMB (per MB)", side = sidePos, line = linePos)
    }
  }
  if (logscale) {
    lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
      segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                        Median_Mutations_log10], y1 = tcga.cohort.med[i, 
                                                                                                      Median_Mutations_log10], col = medianCol)
    })
  }
  else {
    lapply(seq_len(nrow(tcga.cohort.med)), function(i) {
      segments(x0 = i - 1, x1 = i, y0 = tcga.cohort.med[i, 
                                                        Median_Mutations], y1 = tcga.cohort.med[i, Median_Mutations], 
               col = medianCol)
    })
  }
  tcga.cohort = data.table::rbindlist(l = tcga.cohort)
  tcga.cohort[, `:=`(plot_total, NULL)]
  tcga.cohort[, `:=`(TCGA, NULL)]
  pt.test.pval = pt.test.pval[!is.na(Pval)][order(Pval, decreasing = FALSE)]
  list(median_mutation_burden = tcga.cohort.med, mutation_burden_perSample = tcga.cohort, 
       pairwise_t_test = pt.test.pval)
}

