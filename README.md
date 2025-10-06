# MATAGE
MATAGE (additive Mixed-model Association Tests for multi-Allelic Genetic Effects) is a statistical method utilizing a generalized additive mixed model (GAMM) to test for nonlinear genetic effects of Tandem Repeats (TRs), accounting for sample relatedness.

<br />

## Command Line Options
```
Required Arguments:
  --pheno_infile      Path to the phenotype data file.
  --geno_infile       Path to the genotype data file.
  --outfile_dir       Directory to save result files.
  --outfile_name      Basename for output files.
  --ID_name           Column name for subject IDs in the phenotype file.
  --formula           Model formula for the null model (e.g., "Phenotype~Age+Sex"). Must be quoted.
  --family            Model family type: 'linear' or 'logistic'.
  --geno_start_col    Integer specifying the starting column of genotype data in geno_infile.

Optional Arguments:
  --kins_infile       Path to the kinship matrix file. (Default: NULL)
  --total_rows        Total number of variants (rows) to process in geno_infile. (Default: line count of geno_infile, minus header)
  --num_rows          Number of rows to process per job. (Default: same as --total_rows)
  --max_jobs          Maximum number of parallel Rscript jobs to run. (Default: 1)
  -h, --help          Display this help message and exit.
```

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
