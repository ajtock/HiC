#!/applications/R/R-3.3.2/bin/Rscript

# Usage:
# ./HiCdatR_commandArgs.R 10000 10kb 100

args <- commandArgs(trailingOnly = T)
binSize <- 10000
binName <- "10kb"
smth <- 100
binSize <- as.numeric(args[1])
binName <- as.character(args[2])
smth <- as.numeric(args[3])

# Load the necessary functions (general and organism-specific)
library(HiCdatR)
f.source.organism.specific.code("/home/ajt200/HiC/HiCdat/HiCdatR/HiCdat-A-thaliana-TAIR10.R")

# Specify path to directory containing input files and plot output directory
inDir <- "/home/ajt200/HiC/Feng_Jacobsen_2014_Mol_Cell/"
plotDir <- "/home/ajt200/HiC/Feng_Jacobsen_2014_Mol_Cell/plots/"

# Specify the samples and their reduced matrices for R
sampleList <- list(wt = c(paste0("HiC_HindIII_Col0_SRR1504819_read_pairs_mapped_", binName, "_reducedMatrixR.txt")),
                   kss = c(paste0("HiC_HindIII_suvh456_SRR1504827_read_pairs_mapped_", binName, "_reducedMatrixR.txt")))

# Load all samples into a list,
# acccess an individual HiC matrix with binMatList[[sampleName]]
# or binMatList$sampleName
binMatList <- f.load.samples(dataDir = inDir,
                             sampleToFiles = sampleList,
                             binSize = binSize,
                             repetitions = 50)

annotationFromFragments <- f.read.annotation(paste0(inDir, "TAIR10_fragments_", binName, "_ann.txt"),
                                             binSize)

wt_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList[["wt"]],
                                                      binSize = binSize,
                                                      rDir = plotDir,
                                                      outfilePrefix = paste0("HiC_HindIII_Col0_SRR1504819_PCA_", binName),
                                                      annotation = annotationFromFragments)
save(wt_PCA, file = paste0(inDir, "HiC_HindIII_Col0_SRR1504819_PCA_", binName, ".RData"))

kss_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList[["kss"]],
                                                       binSize = binSize,
                                                       rDir = plotDir,
                                                       outfilePrefix = paste0("HiC_HindIII_suvh456_SRR1504827_PCA", binName),
                                                       annotation = annotationFromFragments)
save(kss_PCA, file = paste0(inDir, "HiC_HindIII_suvh456_SRR1504827_PCA_", binName, ".RData"))

load(file = paste0(inDir, "HiC_HindIII_Col0_SRR1504819_PCA_", binName, ".RData"))
load(file = paste0(inDir, "HiC_HindIII_suvh456_SRR1504827_PCA_", binName, ".RData"))

wt_PCA_regToPC_allChrs <- NULL
for(i in 1:5) {
  wt_PCA_regToPC_allChrs <- data.frame(rbind(wt_PCA_regToPC_allChrs, wt_PCA$regToPC[[i]]))
}
kss_PCA_regToPC_allChrs <- NULL
for(i in 1:5) {
  kss_PCA_regToPC_allChrs <- data.frame(rbind(kss_PCA_regToPC_allChrs, kss_PCA$regToPC[[i]]))
}

annotationFromFragments_eigval_wt <- merge(x = annotationFromFragments,
                                           y = wt_PCA_regToPC_allChrs,
                                           by.x = "fragmentNumber",
                                           by.y = "validBins",
                                           all.x = TRUE)
annotationFromFragments_eigval_wt_kss <- merge(x = annotationFromFragments_eigval_wt,
                                               y = kss_PCA_regToPC_allChrs,
                                               by.x = "fragmentNumber",
                                               by.y = "validBins",
                                               all.x = TRUE) 
colnames(annotationFromFragments_eigval_wt_kss) <- c(colnames(annotationFromFragments_eigval_wt_kss)[1:dim(annotationFromFragments)[[2]]],
                                                     "wt_orientedPrincp", "kss_orientedPrincp")

## Plot genome profiles using the Hi-C eigenvalues, which correspond to the
# first principal component of the Pearson's matrix
# The sign of the eigenvalue indicates the compartment:
# (in annotationFromFragments_eigval_wt_kss, A/B [euchromatic/heterochromatic] = +/- ;
#  multiply eigenvalues by -1 so that A/B [euchromatic/heterochromatic] = -/+)
annotationFromFragments_eigval_wt_kss_colnames <- colnames(annotationFromFragments_eigval_wt_kss)
annotationFromFragments_eigval_wt_kss_invertSign <- cbind(annotationFromFragments_eigval_wt_kss,
                                                          annotationFromFragments_eigval_wt_kss$wt_orientedPrincp*-1,
                                                          annotationFromFragments_eigval_wt_kss$kss_orientedPrincp*-1)
colnames(annotationFromFragments_eigval_wt_kss_invertSign) <- c(annotationFromFragments_eigval_wt_kss_colnames,
                                                                "wt_orientedPrincpInvertSign",
                                                                "kss_orientedPrincpInvertSign")

# Genomic definitions
chrs <- c("1", "2", "3", "4", "5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes 
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

annotationFromFragments_eigval_wt_kss_invertSign_cumStart <- NULL
annotationFromFragments_eigval_wt_kss_invertSign_cumEnd <- NULL
for(x in 1:5) {
  annotationFromFragments_eigval_wt_kss_invertSign_chr <- annotationFromFragments_eigval_wt_kss_invertSign[annotationFromFragments_eigval_wt_kss_invertSign$chrom == chrs[x],]
  annotationFromFragments_eigval_wt_kss_invertSign_chr_cumStart <- annotationFromFragments_eigval_wt_kss_invertSign_chr$start+sumchr[x]+1
  annotationFromFragments_eigval_wt_kss_invertSign_cumStart <- c(annotationFromFragments_eigval_wt_kss_invertSign_cumStart,
                                                                 annotationFromFragments_eigval_wt_kss_invertSign_chr_cumStart)
  annotationFromFragments_eigval_wt_kss_invertSign_chr_cumEnd <- annotationFromFragments_eigval_wt_kss_invertSign_chr$end+sumchr[x]+1
  annotationFromFragments_eigval_wt_kss_invertSign_cumEnd <- c(annotationFromFragments_eigval_wt_kss_invertSign_cumEnd,
                                                               annotationFromFragments_eigval_wt_kss_invertSign_chr_cumEnd)
}

# smooth coverage values with MA filter
ma <- rep(1, smth)/smth

filt_vals <- NULL
#filt_vals_noNA <- NULL
for(x in 1:5) {
  df_chr <- annotationFromFragments_eigval_wt_kss_invertSign[annotationFromFragments_eigval_wt_kss_invertSign$chrom == chrs[x],] 
  filt_vals_chr <- NULL
#  filt_vals_chr_noNA <- NULL
  for(y in 5:dim(df_chr)[[2]]) {
     filt_vals_chr_column <- stats::filter(df_chr[,y], filter = ma)
     which_na <- which(is.na(filt_vals_chr_column) == TRUE)
     left_na <- which_na[which(which_na < 100)] 
     left_val <- filt_vals_chr_column[left_na[length(left_na)]+1]
     filt_vals_chr_column[left_na] <- left_val
     right_na <- which_na[which(which_na > 100)]
     right_val <- filt_vals_chr_column[right_na[1]-1]
     filt_vals_chr_column[right_na] <- right_val
#     filt_vals_chr_column_noNA <- filt_vals_chr_column[!is.na(filt_vals_chr_column)]
     filt_vals_chr <- cbind(filt_vals_chr, filt_vals_chr_column)
#     filt_vals_chr_noNA <- cbind(filt_vals_chr_noNA, filt_vals_chr_column_noNA)
  } 
  filt_vals <- rbind(filt_vals, filt_vals_chr)
#  filt_vals_noNA <- rbind(filt_vals_noNA, filt_vals_chr_noNA)
}

annotationFromFragments_eigval_wt_kss_invertSign_colnames <- colnames(annotationFromFragments_eigval_wt_kss_invertSign)
annotationFromFragments_eigval_wt_kss_invertSign_cumCoords_filt <- cbind(annotationFromFragments_eigval_wt_kss_invertSign,
                                                                         annotationFromFragments_eigval_wt_kss_invertSign_cumStart,
                                                                         annotationFromFragments_eigval_wt_kss_invertSign_cumEnd,
                                                                         filt_vals)
colnames(annotationFromFragments_eigval_wt_kss_invertSign_cumCoords_filt) <- c(annotationFromFragments_eigval_wt_kss_invertSign_colnames,
                                                                               "cumStart", 
                                                                               "cumEnd",
                                                                               "den_wt_REC8_HA_Rep1_filt",
                                                                               "den_kss_REC8_HA_Rep1_filt",
                                                                               "wt_orientedPrincp_filt",
                                                                               "kss_orientedPrincp_filt",
                                                                               "wt_orientedPrincpInvertSign_filt",
                                                                               "kss_orientedPrincpInvertSign_filt")

df <- annotationFromFragments_eigval_wt_kss_invertSign_cumCoords_filt

# Function to plot genome-scale coverage of one dataset overlaid with other dataset 
oneVtwoGenomePlot <- function(xplot, one_norm, two_norm, three_norm, four_norm, oneYlabel, twoYlabel) {
  plot(xplot, two_norm, type = "l", lwd = 1.5, col = "darkgreen",
       ylim = c(min(c(two_norm, four_norm), na.rm = T),
                max(c(two_norm, four_norm), na.rm = T)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(one_norm, two_norm, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  axis(side = 4, at = pretty(c(two_norm, four_norm)), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.2, adj = c(0.5, -3.5), labels = twoYlabel, xpd = NA, srt = -90, col = "darkgreen")
#  mtext(side = 4, line = 2.5, cex = 1, text = twoYlabel, col = "darkgreen")
  par(new = T)
  plot(xplot, one_norm, type = "l", lwd = 1.5, col = "magenta",
       ylim = c(min(c(one_norm, three_norm), na.rm = T),
                max(c(one_norm, three_norm), na.rm = T)),
       xlab = "",
       ylab = "")
  axis(side = 2, at = pretty(c(one_norm, three_norm)), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = oneYlabel, col = "magenta")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

# Function to plot genome-scale coverage of one dataset overlaid with other dataset 
oneVtwoGenomePlot2 <- function(xplot, one_norm, two_norm, three_norm, four_norm, oneYlabel, twoYlabel) {
  plot(xplot, one_norm, type = "l", lwd = 1.5, col = "magenta",
       ylim = c(min(c(one_norm, three_norm), na.rm = T),
                max(c(one_norm, three_norm), na.rm = T)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(one_norm, two_norm, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  axis(side = 4, at = pretty(c(one_norm, three_norm)), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1.2, adj = c(0.5, -3.5), labels = oneYlabel, xpd = NA, srt = -90, col = "magenta")
#  mtext(side = 4, line = 2.5, cex = 1, text = twoYlabel, col = "magenta")
  par(new = T)
  plot(xplot, two_norm, type = "l", lwd = 1.5, col = "darkgreen",
       ylim = c(min(c(two_norm, four_norm), na.rm = T),
                max(c(two_norm, four_norm), na.rm = T)),
       xlab = "",
       ylab = "")
  axis(side = 2, at = pretty(c(two_norm, four_norm)), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = twoYlabel, col = "darkgreen")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

# Function to plot genome-scale profile of wt vs mutant (one Y-axis) 
wtVmutantGenomePlot <- function(xplot,
                                wt,
                                mutant,
                                Ylabel,
                                legendLoc,
                                legendLabs,
                                wtCol,
                                mutantCol) {
  plot(xplot, wt, type = "l", lwd = 1.5, col = wtCol,
       ylim = c(min(c(wt, mutant), na.rm = T),
                max(c(wt, mutant), na.rm = T)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(wt, mutant, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, mutant, type = "l", lwd = 1.5, col = mutantCol)
  axis(side = 2, at = pretty(c(wt, mutant)), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(wtCol, mutantCol),
         text.col = c(wtCol, mutantCol),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

pdf(file = paste0(plotDir, "log2_wt_kss_REC8_HA_Rep1_wt_kss_HiC_eigenvalues_genomeplot_", binName, "_rho_v100818.pdf"),
    height = 6.25, width = 18.75)
par(mfrow = c(2, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_wt_REC8_HA_Rep1,
                  two_norm = df$wt_orientedPrincpInvertSign,
                  three_norm = df$den_kss_REC8_HA_Rep1,
                  four_norm = df$kss_orientedPrincpInvertSign,
                  oneYlabel = "wt REC8-HA Rep1",
                  twoYlabel = "wt Hi-C eigenvalue (inverted sign)")
oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_kss_REC8_HA_Rep1,
                  two_norm = df$kss_orientedPrincpInvertSign,
                  three_norm = df$den_wt_REC8_HA_Rep1,
                  four_norm = df$wt_orientedPrincpInvertSign,
                  oneYlabel = "kss REC8-HA Rep1",
                  twoYlabel = "kss Hi-C eigenvalue (inverted sign)")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$wt_orientedPrincpInvertSign,
                    mutant = df$kss_orientedPrincpInvertSign,
                    Ylabel = "Hi-C eigenvalue (inverted sign)",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$den_wt_REC8_HA_Rep1,
                    mutant = df$den_kss_REC8_HA_Rep1,
                    Ylabel = "REC8-HA Rep1",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
dev.off()

pdf(file = paste0(plotDir, "filt", as.character(smth), "_log2_wt_kss_REC8_HA_Rep1_wt_kss_HiC_eigenvalues_genomeplot_", binName, "_rho_v100818.pdf"),
    height = 6.25, width = 18.75)
par(mfrow = c(2, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_wt_REC8_HA_Rep1_filt,
                  two_norm = df$wt_orientedPrincpInvertSign,
                  three_norm = df$den_kss_REC8_HA_Rep1_filt,
                  four_norm = df$kss_orientedPrincpInvertSign,
                  oneYlabel = "wt REC8-HA Rep1",
                  twoYlabel = "wt Hi-C eigenvalue (inverted sign)")
oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_kss_REC8_HA_Rep1_filt,
                  two_norm = df$kss_orientedPrincpInvertSign,
                  three_norm = df$den_wt_REC8_HA_Rep1_filt,
                  four_norm = df$wt_orientedPrincpInvertSign,
                  oneYlabel = "kss REC8-HA Rep1",
                  twoYlabel = "kss Hi-C eigenvalue (inverted sign)")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$wt_orientedPrincpInvertSign,
                    mutant = df$kss_orientedPrincpInvertSign,
                    Ylabel = "Hi-C eigenvalue (inverted sign)",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$den_wt_REC8_HA_Rep1_filt,
                    mutant = df$den_kss_REC8_HA_Rep1_filt,
                    Ylabel = "REC8-HA Rep1",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
dev.off()

pdf(file = paste0(plotDir, "filt", as.character(smth), "_log2_wt_kss_REC8_HA_Rep1_filt", as.character(smth), "_wt_kss_HiC_eigenvalues_genomeplot_", binName, "_rho_v100818.pdf"),
    height = 6.25, width = 18.75)
par(mfrow = c(2, 2))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_wt_REC8_HA_Rep1_filt,
                  two_norm = df$wt_orientedPrincpInvertSign_filt,
                  three_norm = df$den_kss_REC8_HA_Rep1_filt,
                  four_norm = df$kss_orientedPrincpInvertSign_filt,
                  oneYlabel = "wt REC8-HA Rep1",
                  twoYlabel = "wt Hi-C eigenvalue (inverted sign)")
oneVtwoGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                  one_norm = df$den_kss_REC8_HA_Rep1_filt,
                  two_norm = df$kss_orientedPrincpInvertSign_filt,
                  three_norm = df$den_wt_REC8_HA_Rep1_filt,
                  four_norm = df$wt_orientedPrincpInvertSign_filt,
                  oneYlabel = "kss REC8-HA Rep1",
                  twoYlabel = "kss Hi-C eigenvalue (inverted sign)")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$wt_orientedPrincpInvertSign_filt,
                    mutant = df$kss_orientedPrincpInvertSign_filt,
                    Ylabel = "Hi-C eigenvalue (inverted sign)",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
wtVmutantGenomePlot(xplot = df$cumStart+((df$cumEnd-df$cumStart)/2),
                    wt = df$den_wt_REC8_HA_Rep1_filt,
                    mutant = df$den_kss_REC8_HA_Rep1_filt,
                    Ylabel = "REC8-HA Rep1",
                    legendLoc = "topleft",
                    legendLabs = c("Wild type", expression(italic("kyp suvh5 suvh6"))),
                    wtCol = "blue",
                    mutantCol = "red")
dev.off()

