# MATAGE
MATAGE (additive Mixed-model Association Tests for multi-Allelic Genetic Effects) is a statistical method utilizing a generalized additive mixed model (GAMM) to test for nonlinear genetic effects of Tandem Repeats (TRs), accounting for sample relatedness.

<br />

## Command Line Options

<details>
     <summary> <b>List of Options</b> </summary>

```
   --work_dir         
     A working directory.

   --log_dir
     Directory for storing the log files.

   --pheno_infile
     Directory containing the phenotype files.

   --geno_infile
     Directory containing the genotype files.

   --geno_start_col
     The column index at which the genotype data begins in the genotype file.

   --kins_infile
     Directory containing the kinship files. If no kinship matrix is required, set this to NULL.

   --outfile_dir
     Directory for storing the test result output files (not including the file name).

   --name
     The output file name, not including the file extension.

   --ID_name
     Name of the sample ID column in the phenotype file.

   --formula
     The model for the null hypothesis, including the phenotype and covariates but excluding genetic effects.

   --family
     	A description of the distribution and link function to be used in the model. Either "gaussian(link = "identity")" for continuous outcome or "binomial(link = "logit")" for binary outcome.

   --null_model_file
     Directory for storing the null model RData files.

   --total_rows
     Total number of rows (genotypes) in the genotype file.

   --num_rows
     Number of rows (genotypes) to test per job.

   --max_jobs
     Maximum number of jobs to submit simultaneously.

```
</details>
