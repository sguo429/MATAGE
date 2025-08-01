library(data.table)
library(Matrix)
library(GMMAT)

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
pheno_infile <- args[1]
kins_infile <- args[2]
ID_name <- args[3]
formula <- args[4]
null_model_file <- args[5]

# Read in the phenotype data
pheno_data <- fread(file = pheno_infile, header = TRUE, data.table = FALSE)

# Create kinship matrix if exist
if(kins_infile != "NULL"){
    ped <- read.table(kins_infile, header = TRUE)
    N <- nrow(pheno_data)
    idx1 <- match(ped$ID1, pheno_data[ , ID_name])
    idx2 <- match(ped$ID2, pheno_data[ , ID_name])
    idxmin <- pmin(idx1, idx2)
    idxmax <- pmax(idx1, idx2)
    kins <- 0.5*Diagonal(N)+sparseMatrix(i = idxmin, j = idxmax, x = ped$Kinship, dims = c(N,N), symmetric = TRUE)
    colnames(kins) <- rownames(kins) <- pheno_data[ , ID_name]
} else {
    kins <- NULL
}

# Fit the null model
null_model <- glmmkin(
	fixed = formula,   
	data = pheno_data,
	kins = kins,
	id = ID_name,
	family = gaussian(link = "identity"))

save(null_model, file = null_model_file)
