# lieber_project

Usage: call the scripts in the following order:
# read the subjects table in csv format and save it as csv_data.mat
read_csv.m
# read the raw snp data (also load csv_data.mat), and save it as raw_data.mat
read_raw.m
# read the fmri data (also load csv_data.mat), and save it as fmri_data.mat
fetch_fmri.m
# run the main part of the experiment
pipeline_test.m
