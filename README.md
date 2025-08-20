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
     Directory for the phenotype file.

   --geno_infile
     Directory for the genotype file.

   --geno_start_col
     The column index at which the genotype data begins in the genotype file.

   --kins_infile
     Directory containing the kinship file. If no kinship matrix is required, set this to NULL.

   --outfile_dir
     Directory for storing the test result output files (not including the file name).

   --name
     The output file name, not including the file extension.

   --ID_name
     Name of the sample ID column in the phenotype file.

   --formula
     The model for the null hypothesis, including the phenotype and covariates but excluding genotype.

   --family
     A description of the distribution and link function to be used in the model. Either "gaussian(link = "identity")" for continuous outcome or "binomial(link = "logit")" for binary outcome.

   --null_model_file
     Directory for the null model RData file.

   --total_rows
     Total number of rows (genotypes) in the genotype file.

   --num_rows
     Number of rows (genotypes) to test per job.

   --max_jobs
     Maximum number of jobs to submit simultaneously.

```
</details>

<br /> 

## Input Files

MATAGE accepts plain text genotype file, either original ".txt" or compressed ".txt.gz" file.

<br />

## Output File

MATAGE will write results to the output file.
Below are details of the column headers in the output file.

```diff 
ID              - The genotype ID.
n               - The sample size.
knots           - The unique genotype counts.
mean_G          - The mean value of the genotype.
var_G           - The variance of the genotype.
p_value_l       - The linear p-value.
p_value_nl      - The nonlinear p-value.
p_value_joint   - The MATAGE p-value.
```

<br />
