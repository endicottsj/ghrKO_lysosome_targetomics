# Load necessary tools
library(Biobase)
library(limma)
library(readxl)
library(ggplot2)
library(naniar)
library(tidyverse)
library(imputeLCMD)
library(PCAtools)

# Load necessary data
WLP_E <- read_excel("~/Documents/Science/Proteomics Spread Sheets/WLP_Removed_1HW.xlsx", sheet = "Expression Data", col_names = TRUE, col_types = "numeric")
WLP_P <- read_excel("~/Documents/Science/Proteomics Spread Sheets/WLP_Removed_1HW.xlsx", sheet = "Phenotype Data", col_names = TRUE, col_types = "text")
WLP_F <- read_excel("~/Documents/Science/Proteomics Spread Sheets/WLP_Removed_1HW.xlsx", sheet = "Feature Data", col_names = TRUE, col_types = "text")
SN <- unlist(read_excel("~/Documents/Science/Proteomics Spread Sheets/WLP_Removed_1HW.xlsx", sheet = "Sample Names", col_names = FALSE, col_types = NULL))
FN <- unlist(read_excel("~/Documents/Science/Proteomics Spread Sheets/WLP_Removed_1HW.xlsx", sheet = "Feature Names", col_names = FALSE, col_types = NULL))


# replace values of "0" with NA
WLP_E_clean <- replace_with_na_all(WLP_E, ~.x == 0)


# visualize missing values before data imputation
vis_miss(WLP_E_clean[,1:24]) +
  ylab("Proteins") +
  ggtitle("Whole liver proteomics: missing values") +
  theme(plot.margin = margin(10,80,10,10),
        text = element_text(size = 16))
vis_miss(WLP_E_clean[,1:24], cluster = TRUE) +
  ylab("Proteins") +
  ggtitle("Whole liver proteomics: missing values, clustered") +
  theme(plot.margin = margin(10,80,10,10),
        text = element_text(size = 16))


# Determine if distribution of missingness is related to detection counts
WLP_E_shadow <- bind_shadow(WLP_E_clean)

ggplot(WLP_E_shadow, aes(x = rowMeans(WLP_E_clean[,1:24]), color = F.WT.1_NA)) +
  scale_color_manual(values = c("#512DA8", "#E53935"),
                     name = "Legend",
                     labels = c("Not Missing", "Missing")) +
  geom_density() +
  scale_x_log10(limits = c(1e1, 1e12)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  xlab("Average Read Intensity") +
  ylab("Density") +
  theme(axis.line.x.bottom = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        legend.key = element_rect(fill = "white"),
        text = element_text(size = 16),
        panel.background = element_rect(fill = "white", color = NULL))
# Finding that "missingness" is linked to low detection counts, 
# we can use left-censored missing data (LCMD) imputation, after nomalization
# Variance Stabilizing Normalization works but is not the accepted standard in the field
# WLP_E_VSN <- normalizeVSN(WLP_E_clean)


# LIMMA groundwork: create Expression Set Object
mat_WLP_E <- data.matrix(WLP_E_clean, rownames.force = FALSE)
rownames(mat_WLP_E) <- FN
rownames(WLP_P) <- SN
rownames(WLP_F) <- FN

WLPset <- ExpressionSet(assayData = mat_WLP_E,
                        phenoData = AnnotatedDataFrame(WLP_P),
                        featureData = AnnotatedDataFrame(WLP_F))

#Quantile Normalization: check initial distribution (optional)
exprs(WLPset) <-log(exprs(WLPset))
plotDensities(WLPset, legend = FALSE)
WLPset_norm <- WLPset
exprs(WLPset_norm) <- normalizeBetweenArrays(exprs(WLPset_norm))
plotDensities(WLPset_norm, legend = FALSE)

#Optional Print Data
print_WLPset_norm <- as.data.frame(exprs(WLPset_norm), row.names = FN, col.names = SN)
write.csv(print_WLPset_norm, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_QN.csv")


# Data imputation:
WLPset_QN_I <- WLPset_norm
exprs(WLPset_QN_I) <- impute.MinProb(exprs(WLPset_QN_I))


# Data imputation quality control, do not need to run normally:
WLP_imputed <- impute.MinProb(exprs(WLPset_QN_I))
WLP_imputed_df <- as.data.frame(WLP_imputed, row.names = FN, col.names = SN)
write.csv(WLP_imputed_df, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_imputed_df.csv")
vis_miss(WLP_imputed_df[,1:24]) +
  ylab("Proteins") +
  ggtitle("Whole liver proteomics: post-imputation") +
  theme(plot.margin = margin(10,80,10,10),
        text = element_text(size = 16))


#  LIMMA groundwork: name a linear model
study_design <- model.matrix(~0 + Group, data = pData(WLPset_QN_I))
fit_1 <- lmFit(WLPset_QN_I, study_design)


# WLP Female Genotype Comparison Gx v WT
F_c_matrix <- makeContrasts(F_GxvWT = GroupF.Gx - GroupF.WT,
                            levels = study_design)
fit_FGC <- contrasts.fit(fit_1, contrasts = F_c_matrix)
fit_FGC_eBayes <- eBayes(fit_FGC)
# check that it worked
volcanoplot(fit_FGC_eBayes, highlight = 10)
WLP_FGC <- topTable(fit_FGC_eBayes, number = nrow(fit_FGC_eBayes), sort.by = "none")
write.csv(WLP_FGC, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_F_geno_contrast.csv")


# WLP Male Genotype Comparison Gx v WT
M_c_matrix <- makeContrasts(M_GxvWT = GroupM.Gx - GroupM.WT,
                            levels = study_design)
fit_MGC <- contrasts.fit(fit_1, contrasts = M_c_matrix)
fit_MGC_eBayes <- eBayes(fit_MGC)
# check that it worked
volcanoplot(fit_MGC_eBayes, highlight = 10)
WLP_MGC <- topTable(fit_MGC_eBayes, number = nrow(fit_MGC_eBayes), sort.by = "none")
write.csv(WLP_MGC, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_M_geno_contrast.csv")


# WLP Pooled, Genotype Comparison Gx v WT, Ignoring Sex
pooled_study_design <- model.matrix(~0 + Genotype, data = pData(WLPset_QN_I))
pooled_fit <- lmFit(WLPset_QN_I, pooled_study_design)
FM_c_matrix <- makeContrasts(FM_GxvWT = GenotypeGx - GenotypeWT,
                             levels = pooled_study_design)
fit_FMGC <- contrasts.fit(pooled_fit, contrasts = FM_c_matrix)
fit_FMGC_eBayes <- eBayes(fit_FMGC)
# check that it worked
volcanoplot(fit_FMGC_eBayes, highlight = 10)
WLP_FMGC <- topTable(fit_FMGC_eBayes, number = nrow(fit_FMGC_eBayes), sort.by = "none")
write.csv(WLP_FMGC, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_pooled_geno_contrast.csv")

# F v M comparisons for controls and for Gx
C_s_matrix <- makeContrasts(C_FvM = GroupM.WT - GroupF.WT,
                            levels = study_design)
fit_Cs <- contrasts.fit(fit_1, contrasts = C_s_matrix)
fit_Cs_eBayes <- eBayes(fit_Cs)
# check that it worked
volcanoplot(fit_Cs_eBayes, highlight = 10)
WLP_Cs <- topTable(fit_Cs_eBayes, number = nrow(fit_Cs_eBayes), sort.by = "none")
write.csv(WLP_Cs, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_C_Sex_contrast.csv")

# Now for Gx
G_s_matrix <- makeContrasts(G_FvM = GroupM.Gx - GroupF.Gx,
                            levels = study_design)
fit_Gs <- contrasts.fit(fit_1, contrasts = G_s_matrix)
fit_Gs_eBayes <- eBayes(fit_Gs)
# check that it worked
volcanoplot(fit_Gs_eBayes, highlight = 10)
WLP_Gs <- topTable(fit_Gs_eBayes, number = nrow(fit_Gs_eBayes), sort.by = "none")
write.csv(WLP_Gs, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Gx_Sex_contrast.csv")

# Code for pathway analysis

#Females only
FGC_GO <- goana(fit_FGC_eBayes, geneid = "Entrez", species = "Mm")
write.csv(FGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Female_GO.csv")
# No FDR
FGC_GO <- goana(fit_FGC_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(FGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Female_GO_NoFDR.csv")

#Males only
MGC_GO <- goana(fit_MGC_eBayes, geneid = "Entrez", species = "Mm")
write.csv(MGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Male_GO.csv")
# No FDR
MGC_GO <- goana(fit_MGC_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(MGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Male_GO_NoFDR.csv")


#Pooled Data
PooledGC_GO <- goana(fit_FMGC_eBayes, geneid = "Entrez", species = "Mm")
write.csv(PooledGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Pooled_GO.csv")
# No FDR
PooledGC_GO <- goana(fit_FMGC_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(PooledGC_GO, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_Pooled_GO_NoFDR.csv")

# Two factor design
contrast_matrix <- makeContrasts(Interaction = (GroupM.Gx - GroupF.Gx) - (GroupM.WT - GroupF.WT),
                                 levels = study_design)

fit_Int <- contrasts.fit(fit_1, contrasts = contrast_matrix)
fit_Int_eBayes <- eBayes(fit_Int)
results <- decideTests(fit_Int_eBayes)
summary(results)
# Getting the Data
Two_Factor_WLP <- topTable(fit_Int_eBayes, number = nrow(fit_Int_eBayes), sort.by = "none")
volcanoplot(fit_Int_eBayes, highlight = 10)
write.csv(Two_Factor_WLP, file = "~/Documents/Science/Proteomics Spread Sheets/WLP_2F.csv")
