#!/bin/bash
ml Miniconda3
source activate pyclone #point to the Pyclone directory path
PyClone build_mutations_file --in_file sample_data.txt --out_file sample_data.yaml --prior total_copy_number
PyClone run_analysis --config_file sample_config.yaml
PyClone build_table --burnin 1000 --config_file sample_config.yaml --out_file sample.loci.table.txt --table_type loci
PyClone build_table --burnin 1000 --config_file sample_config.yaml --out_file sample.cluster.table.txt --table_type cluster
source deactivate
