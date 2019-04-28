# Libraries needed
library(limma)

# Import Acetylproteins
input <- read.csv("limma_export_dkotac_proteins.csv", 
                            row.names=1)

# Convert to log2
log2 <- log2(input)


# Design linear model with no intercept and no interaction
group <- factor(c("DKO_Sham", "DFC_Sham", "DFC_TAC", 
                  "DFC_Sham", "DKO_TAC", "DFC_TAC", 
                  "DKO_TAC", "DKO_TAC", "DFC_TAC", 
                  "DKO_Sham", "DFC_Sham"))

design <- model.matrix(~0+group)

# Rename column names in the design
colnames(design) <- c("DFC_Sham", "DFC_TAC", "DKO_Sham", "DKO_TAC")

# Limma Step 1: Least Squares Estimates 
fit <- lmFit(log2, design)

# Generate the contrast map
contMat <- makeContrasts(DKO_TACvsDFC_TAC = DKO_TAC-DFC_TAC,
                         DKO_ShamvsDFC_Sham = DKO_Sham-DFC_Sham,
                         DKO_ShamvsDFC_TAC = DKO_Sham-DFC_TAC,
                         DFC_TACvsDFC_Sham = DFC_TAC-DFC_Sham,
                         DKO_TACvsDKO_Sham = DKO_TAC-DKO_Sham,
                         level = design)

# Fit the Contrasts to the product of lmFit
fit.cont <- contrasts.fit(fit, contMat)

# Now perform the eBayes on this new fit
fit.cont.eb <- eBayes(fit.cont)


# Perform the comparisons, with each contrast individually
# Use Benj. Hoch. for With FDR = 0.05
fit.cont.eb.decide <- decideTests(fit.cont.eb, 
                                  method="separate", 
                                  adjust.method = "BH", 
                                  p.value = 0.05)

# Export the data
write.fit(fit.cont.eb, 
          fit.cont.eb.decide, 
          "eb_fit_dkotac_proteins.tsv")