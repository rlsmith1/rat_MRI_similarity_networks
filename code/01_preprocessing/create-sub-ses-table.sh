#!/bin/bash

base_dir=/data/CamRat/WHS_preprocessing
cut --complement -f 3,4 -d, ${base_dir}/scan_table_for_registration.txt | uniq > ${base_dir}/sub_ses.txt 
