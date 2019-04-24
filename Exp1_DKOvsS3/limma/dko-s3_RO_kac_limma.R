# Libraries needed
library(limma)

# Import Acetylpeptides
log2 <- read.csv("limma_export_dko-s3_RO_kac.csv", 
                      row.names=1)

# Don't convert to log2 - it is already!

# Design linear model with no intercept and no interaction
# Group is the Genotype
group <- factor(c('S3KO','S3KO','S3KO',
                  'S3FL','S3FL',
                  'DFC','DFC',
                  'DKO','DKO','DKO'))

design <- model.matrix(~0+group)

# Rename column names in the design
colnames(design) <- c('DFC','DKO', 'S3FL','S3KO')

# Limma Step 1: Least Squares Estimates 
fit <- lmFit(log2, design)

# Generate the contrast map
contMat <- makeContrasts(DKOvsDFC = DKO - DFC,
                         S3KOvsS3FL = S3KO - S3FL,
                         DKOvsS3KO = DKO - S3KO,
                         DFCvsS3FL = DFC - S3FL,
                         levels=design)

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
          "eb_fit_dko-s3_RO_kac.tsv")