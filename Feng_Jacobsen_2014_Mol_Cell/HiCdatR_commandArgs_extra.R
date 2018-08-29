#!/applications/R/R-3.3.2/bin/Rscript

# Usage:
# ./HiCdatR_commandArgs.R 20000 20kb

args <- commandArgs(trailingOnly = T)
#binSize <- 20000
#binName <- "20kb"
binSize <- as.numeric(args[1])
binName <- as.character(args[2])

# Load the necessary functions (general and organism-specific)
library(HiCdatR)
f.source.organism.specific.code("/home/ajt200/HiC/HiCdat/HiCdatR/HiCdat-A-thaliana-TAIR10.R")

# Specify path to directory containing input files and plot output directory
inDir <- "/home/ajt200/HiC/Feng_Jacobsen_2014_Mol_Cell/"
plotDir <- "/home/ajt200/HiC/Feng_Jacobsen_2014_Mol_Cell/plots/"

# Specify the samples and their reduced matrices for R
sampleList <- list(Col = c(paste0("HiC_HindIII_Col0_SRR1504819_read_pairs_mapped_", binName, "_reducedMatrixR.txt")),
                   kss = c(paste0("HiC_HindIII_suvh456_SRR1504827_read_pairs_mapped_", binName, "_reducedMatrixR.txt")))

# Load all samples into a list,
# acccess an individual HiC matrix with binMatList[[sampleName]]
# or binMatList$sampleName
binMatList <- f.load.samples(dataDir = inDir,
                             sampleToFiles = sampleList,
                             binSize = binSize,
                             repetitions = 50)

# Draw a sample correlation matrix
f.HiC.correlation.matrix(dataMatrixList = binMatList,
                         rDir = plotDir,
                         outfile = "HiC_HindIII_Col0_SRR1504819_suvh456_SRR1504827_correlation_matrix",
                         corMethod = "pearson",
                         summaryFunction = median,
                         useOnlyHighVar = TRUE)
#removed 5958 values due to low variance
#removed 6853 values due to low variance
#removed 5958 values due to low variance

# Failed attempt to prevent use of tiff() function (no longer supported in R)
GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT <- TRUE

# Draw heatmaps of the raw intereaction frequencies
# and the correlated interaction frequencies
f.plot.XY.matrix(matrixToPlot = binMatList[["Col"]],
                 binSize = binSize,
                 axStep = 1e5,
                 rDir = plotDir,
                 outfile = "HiC_HindIII_Col0_SRR1504819_raw_interaction_matrix")
#Error in .External2(C_X11, paste("tiff::", comp, ":", filename, sep = ""),  : 
#  unable to start device TIFF
#In addition: Warning message:
#In tiff(file.path(rDir, paste(outfile, ".tiff", sep = "")), width = 2400,  :
#  no tiff support in this version of R

f.plot.XY.matrix(matrixToPlot = binMatList[["Col"]],
                 binSize = binSize,
                 axStep = 1e6,
                 rDir = plotDir,
                 outfile = "HiC_HindIII_Col0_SRR1504819_cor_interaction_matrix",
                 doCor = TRUE)
#Error in .External2(C_X11, paste("tiff::", comp, ":", filename, sep = ""),  : 
#  unable to start device TIFF
#In addition: Warning message:
#In tiff(file.path(rDir, paste(outfile, ".tiff", sep = "")), width = 2400,  :
#  no tiff support in this version of R

# Calculate the differences between two Hi-C samples of interest
## Relative differences
# Differences between two Hi-C samples A and B can be visualised based on the 
# relative difference of interaction frequencies.
# For each matrix entry (i.e., a pixel at row i and column j), the difference
# between the two samples is divided by the average value:
# Rij = (Aij - Bij)/((Aij + Bij)/2)
# A pair of Hi-C samples can be visualised using the function f.plot.relative.difference() .
# Multiple samples are compared to each other with f.compare.samples.relative.difference() .
###f.plot.relative.difference(dataMatrixA = binMatList[["Col"]],
###                           dataMatrixB = binMatList[["kss"]]

# Signed differences (SDM)
###f.plot.signed.difference(binMatList[["Col"]],

annotationFromFragments <- f.read.annotation(paste0(inDir, "TAIR10_fragments_", binName, "_ann.txt"),
                                             binSize)

wt_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList[["Col"]],
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
test <- seq(1, 1000, by = 1)
j = 50
ma <- rep(1, test[j])/test[j]

filt_vals <- NULL
filt_vals_noNA <- NULL
for(x in 1:5) {
  df_chr <- annotationFromFragments_eigval_wt_kss_invertSign[annotationFromFragments_eigval_wt_kss_invertSign$chrom == chrs[x],] 
  filt_vals_chr <- NULL
  filt_vals_chr_noNA <- NULL
  for(y in 5:dim(df_chr)[[2]]) {
     filt_vals_chr_column <- stats::filter(df_chr[,y], filter = ma)
     which_na <- which(is.na(filt_vals_chr_column) == TRUE)
     left_na <- which_na[which(which_na < 100)] 
     left_val <- filt_vals_chr_column[left_na[length(left_na)]+1]
     filt_vals_chr_column[left_na] <- left_val
     right_na <- which_na[which(which_na > 100)]
     right_val <- filt_vals_chr_column[right_na[1]-1]
     filt_vals_chr_column[right_na] <- right_val
     filt_vals_chr_column_noNA <- filt_vals_chr_column[!is.na(filt_vals_chr_column)]
     filt_vals_chr <- cbind(filt_vals_chr, filt_vals_chr_column)
     filt_vals_chr_noNA <- cbind(filt_vals_chr_noNA, filt_vals_chr_column_noNA)
  } 
  filt_vals <- rbind(filt_vals, filt_vals_chr)
  filt_vals_noNA <- rbind(filt_vals_noNA, filt_vals_chr_noNA)
}  

annotationFromFragments_eigval_wt_kss_invertSign_colnames <- colnames(annotationFromFragments_eigval_wt_kss_invertSign)
annotationFromFragments_eigval_wt_kss_invertSign_cumCoords <- cbind(annotationFromFragments_eigval_wt_kss_invertSign,
                                                                    annotationFromFragments_eigval_wt_kss_invertSign_cumStart,
                                                                    annotationFromFragments_eigval_wt_kss_invertSign_cumEnd)
colnames(annotationFromFragments_eigval_wt_kss_invertSign_cumCoords) <- c(annotationFromFragments_eigval_wt_kss_invertSign_colnames, "cumStart", "cumEnd")

# Function to plot genome-scale coverage of one dataset overlaid with other dataset 
oneVtwoGenomePlot <- function(xplot, one_norm, two_norm, oneYlabel, twoYlabel) {
  plot(xplot, two_norm, type = "l", lwd = 1.5, col = "blue",
       ylim = c(min(two_norm, na.rm = T), max(two_norm, na.rm = T)),
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(one_norm, two_norm, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  axis(side = 4, at = pretty(two_norm), lwd.tick = 1.5)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -4.0), labels = twoYlabel, xpd = NA, srt = -90, col = "blue")
#  mtext(side = 4, line = 2.5, cex = 1, text = twoYlabel, col = "blue")
  par(new = T)
  plot(xplot, one_norm, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(one_norm, na.rm = T), max(one_norm, na.rm = T)),
       xlab = "",
       ylab = "")
  axis(side = 2, at = pretty(one_norm), lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = oneYlabel, col = "red")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  box(lwd = 1.5)
}

pdf(file = paste0(plotDir, "log2_wt_REC8_HA_Rep1_wt_HiC_eigenvalues_genomeplot_", binName, "_rho_v090818.pdf"),
    height = 3.125, width = 9.375)
par(mfcol = c(1, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

oneVtwoGenomePlot(xplot = annotationFromFragments_eigval_wt_kss_invertSign_cumCoords$cumStart+((annotationFromFragments_eigval_wt_kss_invertSign_cumCoords$cumEnd-annotationFromFragments_eigval_wt_kss_invertSign_cumCoords$cumStart)/2),
                  one_norm = annotationFromFragments_eigval_wt_kss_invertSign_cumCoords$den_wt_REC8_HA_Rep1,
                  two_norm = annotationFromFragments_eigval_wt_kss_invertSign_cumCoords$orientedPrincpInvertSign,
                  oneYlabel = "wt REC8-HA Rep",
                  twoYlabel = "wt Hi-C eigenvalue (inverted sign)")
dev.off()


