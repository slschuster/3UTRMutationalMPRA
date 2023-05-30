#! usr/bin/env/python3
_author__="Pushpa Itagi"
__purpose__="Generate the required Pyclone (assess clonality) input files from the samples provided. This script needs directory path and input files to generate the required inputs for Pyclone"
__date__ = "05252023"

#Load Modules
import os, sys
from pathlib import Path
import pandas as pd


#Read files and directories
path_files = "path_to_files_mutations" #path to input files, one sample per patient ~183 samples for prostate cancer study
path_tumor_fraction = "file_containing_tumorfraction_perSample" #file containing the tumor fraction and sample details
dir = os.listdir(path_files)
out_dir = "output_directory_path"

f = open("run_pyclone_samples.sh","w")
df = pd.read_csv(path_tumor_fraction,sep="\t")
work_dir = "" #path to your working directory
delimiter = "" #specify a delimiter to match the file pattern, if needed

#For loop to process each file and output the config files
for files in dir:
    sample_id = files.split(delimiter)[0] #extract the sample_id
    file_name = path_files + sample_id
    create_dir = out_dir + sample_id

    if not os.path.exists(create_dir):
        # if the  directory is not present# then create it.
        os.makedirs(create_dir) #create a dirctory for each patient
            print("Directory created in create_dir ",create_dir)

    #generate all the required template files for runing pyclone
    copy_input_file = "../" + files + " "+create_dir +"/"
    config_file = out_dir+sample_id+"/"+sample_id +"_config.yaml"
    copy_str1 = "cp "+  work_dir+"template_config.yaml " + config_file #to copy the template config file to a per sample config file
    template_file = out_dir+sample_id+"/"+sample_id +".sh" to copy the template shell file to a per sample shell file
    copy_str2 = "cp "+  work_dir+"template.sh " + template_file
    os.system(copy_input_file)
    os.system(copy_str1)
    os.system(copy_str2)
    print("Files placed in .....",create_dir)
    purity = df[df['sampleName']==sample_id]['cellularity'].tolist()
    #do all the processing for tumor fraction, sample name and write these files into the new directory created for each sample
    sample_purity = ' '.join(map(str, purity))
    working_dir = create_dir +"/"
    config_temp = Path(config_file)
    config_temp.write_text(config_temp.read_text().replace("999-999",sample_id))
    config_temp.write_text(config_temp.read_text().replace("0.99999",str(sample_purity)))
    config_temp.write_text(config_temp.read_text().replace("working_dir: ./",("working_dir: "+working_dir)))
    shell_temp = Path(template_file)
    shell_temp.write_text(shell_temp.read_text().replace("999-999",sample_id))
    shell_temp.write_text(shell_temp.read_text().replace("..../",(working_dir)))
    print("Config files processed .....",create_dir)
    #Write the path containing the shell script to the file which will have all samples information
    f.write( template_file)
    f.write("\n")
