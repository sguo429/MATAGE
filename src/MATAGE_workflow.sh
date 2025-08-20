#!/bin/bash
work_dir=your_directory

log_dir=${work_dir}/log

pheno_infile=${work_dir}/pheno_family/pheno.txt

geno_infile=${work_dir}/geno.txt

geno_start_col=3

kins_infile=${work_dir}/kinship.txt
#kins_infile=NULL

outfile_dir=${work_dir}/results

name=outfile_name

ID_name=ID

formula="Phenotype~Age+Sex"

family="gaussian(link=\"identity\")"
#family="binomial(link=\"logit\")"

null_model_file=${work_dir}/null_model/${name}_null_model.RData

total_rows=1809

num_rows=362

max_jobs=5

# Function to check and wait for jobs
wait_for_jobs() {
    while [ $(jobs | wc -l) -ge $max_jobs ]; do
        sleep 1
    done
}

cd ${work_dir}

# Fit the null model
Rscript ./scripts/Fit_null_model.R $pheno_infile $kins_infile $ID_name $formula $family $null_model_file

# Loop through the start_row values (from 1 to total_rows with steps of num_rows)
for start_row in $(seq 1 $num_rows $total_rows); do
    wait_for_jobs

    end_row=$((start_row + num_rows - 1))
    outfile=${outfile_dir}/${name}_${start_row}_${end_row}.txt

    # Run the R script
    /usr/bin/time -v Rscript ./scripts/MATAGE_general.R $geno_infile $start_row $num_rows $outfile $null_model_file $geno_start_col \
        > ${log_dir}/${name}_output_${start_row}.log 2>&1 \
        2>> ${log_dir}/${name}_time_${start_row}.log &
done

wait
