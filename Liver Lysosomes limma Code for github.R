# Code for LIMMA analysis of liver lysosomes
# Load required tools
library(Biobase)
library(limma)
library(readxl)
library(ggplot2)
library(naniar)
library(tidyverse)
library(imputeLCMD)
library(stats)
library(PCAtools)

# Load the Data

LLP_E <- read_excel("~/Documents/Science/Proteomics Spread Sheets/LLP_E_Set.xlsx", sheet = "Expression Data", col_names = TRUE, col_types = "numeric")
LLP_P <- read_excel("~/Documents/Science/Proteomics Spread Sheets/LLP_E_Set.xlsx", sheet = "Phenotype Data", col_names = TRUE, col_types = "text")
LLP_F <- read_excel("~/Documents/Science/Proteomics Spread Sheets/LLP_E_Set.xlsx", sheet = "Feature Data", col_names = TRUE, col_types = "text")
SN <- unlist(read_excel("~/Documents/Science/Proteomics Spread Sheets/LLP_E_Set.xlsx", sheet = "Sample Names", col_names = FALSE, col_types = NULL))
FN <- unlist(read_excel("~/Documents/Science/Proteomics Spread Sheets/LLP_E_Set.xlsx", sheet = "Feature Names", col_names = FALSE, col_types = NULL))

# replace values of "0" with NA
LLP_E_clean <- replace_with_na_all(LLP_E, ~.x == 0)
write.csv(LLP_E_clean, file = "~/Documents/Science/proteomics spread sheet output/LLP_E_clean.csv")

# visualize missing values before data imputation
LLP_NA_table <- read_excel("~/Documents/Science/proteomics spread sheet output/LLP_E_clean.xlsx", sheet = "LLP_E_clean", col_names = TRUE, col_types = "numeric")
NA_factor <- factor(LLP_NA_table$NA_group, levels = c(0, 1, 2))

ggplot(LLP_NA_table, aes(x = Row_AVG, color = NA_factor)) +
  geom_density() +
  scale_x_log10(limits = c(1e2, 1e11), breaks = c(1e2, 1e4, 1e6, 1e8, 1e10)) +
  scale_color_discrete(name = "Missing Values", labels = c("None", "1 or 2", "> 2")) +
  xlab("Raw counts, averaged across all samples") +
  ylab("Density") +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 5),
        plot.margin = margin(20,20,20,20),
        text = element_text(size = 15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(color = "#ECEFF1"),
        panel.background = element_rect(color = "black", fill = "white"))

vis_miss(LLP_E_clean[,1:24]) +
  ylab("Proteins") +
  ggtitle("Lysosome proteomics: missing values") +
  theme(plot.margin = margin(10,80,10,10),
        text = element_text(size = 16))
vis_miss(LLP_E_clean[,1:24], cluster = TRUE) +
  ylab("Proteins") +
  ggtitle("Lysosome proteomics: missing values, clustered") +
  theme(plot.margin = margin(10,80,10,10),
        text = element_text(size = 16))

# LIMMA groundwork: create Expression Set Object
mat_LLP_E <- data.matrix(LLP_E_clean, rownames.force = FALSE)
rownames(mat_LLP_E) <- FN
rownames(LLP_P) <- SN
rownames(LLP_F) <- FN

LLPset <- ExpressionSet(assayData = mat_LLP_E,
                        phenoData = AnnotatedDataFrame(LLP_P),
                        featureData = AnnotatedDataFrame(LLP_F))

#Quantile Normalization: check initial distribution (optional)
exprs(LLPset) <-log(exprs(LLPset))
plotDensities(LLPset, legend = FALSE)
LLPset_norm <- LLPset
exprs(LLPset_norm) <- normalizeBetweenArrays(exprs(LLPset_norm))
plotDensities(LLPset_norm, legend = FALSE)

# Data imputation:
LLPset_QN_I <- LLPset_norm
exprs(LLPset_QN_I) <- impute.MinProb(exprs(LLPset_QN_I))

#Optional Print Data
print_LLPset_QN_I <- as.data.frame(exprs(LLPset_QN_I), row.names = FN, col.names = SN)
write.csv(print_LLPset_QN_I, file = "~/Documents/Science/Proteomics Spread Sheets/LLP_QN_I.csv")

# LIMMA groundwork: name a linear model
study_design <- model.matrix(~0 + Group, data = pData(LLPset_QN_I))
fit_1 <- lmFit(LLPset_QN_I, study_design)

# Build a two-factor contrast matrix
contrast_matrix <- makeContrasts(Interaction = (GroupG.Y - GroupG.N) - (GroupC.Y - GroupC.N),
                                 levels = study_design)

fit_Int <- contrasts.fit(fit_1, contrasts = contrast_matrix)
fit_Int_eBayes <- eBayes(fit_Int)
results <- decideTests(fit_Int_eBayes)
summary(results)
# Getting the Data
volcanoplot(fit_Int_eBayes, highlight = 10)
Two_Factor_LLP <- topTable(fit_Int_eBayes, number = nrow(fit_Int_eBayes), sort.by = "none")
write.csv(Two_Factor_LLP, file = "~/Documents/Science/proteomics spread sheet output/LLP_2F.csv")

# Code for pathway analysis
# No FDR
LLP_2F_NoFDR_KEGG <- kegga(fit_Int_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(LLP_2F_NoFDR_KEGG, file = "~/Documents/Science/Proteomics Spread Sheets/LLP_2F_KEGG_NoFDR.csv")
LLP_2F_NoFDR_GO <- goana(fit_Int_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(LLP_2F_NoFDR_GO, file = "~/Documents/Science/Proteomics Spread Sheets/LLP_2F_GO_NoFDR.csv")



# Make a contrast matrix for Genotype effect
study_design2 <- model.matrix(~0 + Genotype, data = pData(LLPset_QN_I))
fit_2 <- lmFit(LLPset_QN_I, study_design2)
Geno_contrast_matrix <- makeContrasts(Genotype_Effect = GenotypeGx - GenotypeWT,
                                 levels = study_design2)
fit_Geno <- contrasts.fit(fit_2, contrasts = Geno_contrast_matrix)
fit_Geno_eBayes <- eBayes(fit_Geno)
# Check that it worked
volcanoplot(fit_Geno_eBayes, highlight = 10)
# Get the data
Geno_LLP <- topTable(fit_Geno_eBayes, number = nrow(fit_PBS_eBayes), sort.by = "none")
write.csv(Geno_LLP, file = "~/Documents/Science/proteomics spread sheet output/Geno_LLP.csv")
# No FDR pathwat analysis
Geno_LLP_NoFDR_KEGG <- kegga(fit_Geno_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(Geno_LLP_NoFDR_KEGG, file = "~/Documents/Science/Proteomics Spread Sheets/Geno_LLP_KEGG_NoFDR.csv")
Geno_LLP_NoFDR_GO <- goana(fit_Geno_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(Geno_LLP_NoFDR_GO, file = "~/Documents/Science/Proteomics Spread Sheets/Geno_LLP_GO_NoFDR.csv")


# Make a contrast matrix for (pooled) Leupeptin Effect
study_design3 <- model.matrix(~0 + Injection, data = pData(LLPset_QN_I))
fit_3 <- lmFit(LLPset_QN_I, study_design3)
Inj_contrast_matrix <- makeContrasts(Injection_Effect = InjectionLeupeptin - InjectionPBS,
                                      levels = study_design3)
fit_Inj <- contrasts.fit(fit_3, contrasts = Inj_contrast_matrix)
fit_Inj_eBayes <- eBayes(fit_Inj)
# Check that it worked
volcanoplot(fit_Inj_eBayes, highlight = 10)
# Get the data
Inj_LLP <- topTable(fit_Inj_eBayes, number = nrow(fit_Inj_eBayes), sort.by = "none")
write.csv(Inj_LLP, file = "~/Documents/Science/proteomics spread sheet output/Inj_LLP.csv")
# No FDR pathwat analysis
Inj_LLP_NoFDR_KEGG <- kegga(fit_Inj_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(Inj_LLP_NoFDR_KEGG, file = "~/Documents/Science/Proteomics Spread Sheets/Inj_LLP_KEGG_NoFDR.csv")
Inj_LLP_NoFDR_GO <- goana(fit_Inj_eBayes, geneid = "Entrez", species = "Mm", FDR = 1)
write.csv(Inj_LLP_NoFDR_GO, file = "~/Documents/Science/Proteomics Spread Sheets/Inj_LLP_GO_NoFDR.csv")



# We need to consider the effects of leupeptin on each genotype separately,
# as a quality control measurement
# i.e. did leupeptin reduce cathepsins in both genotypes? increase sqstm1?
# we'll re-use the original study design from above
study_design <- model.matrix(~0 + Group, data = pData(LLPset_QN_I))
fit_1 <- lmFit(LLPset_QN_I, study_design)

# Build a two-factor contrast matrix
CL_contrast_matrix <- makeContrasts(C_leu = GroupC.Y - GroupC.N,
                                 levels = study_design)
GL_contrast_matrix <- makeContrasts(G_leu = GroupG.Y - GroupG.N,
                                    levels = study_design)

fit_CL <- contrasts.fit(fit_1, contrasts = CL_contrast_matrix)
fit_CL_eBayes <- eBayes(fit_CL)

fit_GL <- contrasts.fit(fit_1, contrasts = GL_contrast_matrix)
fit_GL_eBayes <- eBayes(fit_GL)

# Getting the Data
volcanoplot(fit_CL_eBayes, highlight = 10)
volcanoplot(fit_GL_eBayes, highlight = 10)

CL_Effect_LLP <- topTable(fit_CL_eBayes, number = nrow(fit_CL_eBayes), sort.by = "none")
write.csv(CL_Effect_LLP, file = "~/Documents/Science/proteomics spread sheet output/LLP_CL.csv")

GL_Effect_LLP <- topTable(fit_GL_eBayes, number = nrow(fit_GL_eBayes), sort.by = "none")
write.csv(GL_Effect_LLP, file = "~/Documents/Science/proteomics spread sheet output/LLP_GL.csv")




