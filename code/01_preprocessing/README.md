##########################################################################

## Code to preprocess rat structural MRI data using AFNI software tools ##
## Smith et al Net Neuro 2026					                    	##

##########################################################################

--------------------------------------------------------------------------

00.make-scan-list.sh

 - Generates table of scans and session dates for use in pipeline
 - This is specific to data format on the Cambridge cluster; unlikely to be useful in other contexts
 - Specifications: set out_dir; this is the output path where the table will be generated
 - can run directly in terminal using sh

--------------------------------------------------------------------------

01.data-to-bids.R & run-01.sh

 - Copies data over from original directory to desired output directory in BIDS format <https://bids.neuroimaging.io>
 - run-01.sh will run the R script
 - Specifications: Need to set path in R file & sh file
 - also generates scan_table_for_registration.txt in out_dir, which is a table of all the scans to run through the registration pipeline (including their init_scale)

--------------------------------------------------------------------------

02.align-to-atlas.R & run-02.sh

 - orients scans to align with specified atlas using AFNI tools
 - Specifications: need to update path in both the R and shell script to point to your desired home directory
 - note: make sure AFNI binaries are updates
 - other note: check scan orientation before running 03! The orientation specified is specific to Cambridge scanner output

--------------------------------------------------------------------------

03.align-to-atlas.R & run-03.sh

 - registration script, calls @animal_warper AFNI function
 - note init_scale: this is the parameter to shift after QC runs

--------------------------------------------------------------------------

copy-qc-images.sh

 - Copies QC images (axial, coronal, and sagittal planes for each scan) from derived directory to registration QC folder in the outputs directory
 - Transfer the QC images directory to local space to check registration success
 - IMPORTANT: need to look at each QC image for each scan to check registration success!! If the scan failed, try changing inti_scale. Currently, I create a copy of scan_table_for_registration.txt, save it as an xlsx file called qc{qc_round}_review.xlsx, and use this to record the results from the registration of that scan using that init_scale parameter. I add three columns: 'registration_success' (1 indicates success, 0 indicates failed), 'init_scale2' (suggested inti_scale parameter for next round), and 'notes' (describes what the problem was in the failed scan). I then use this file to create the next scan_table_for_registration.txt, export that back to the cluster and run the next round of registration 

--------------------------------------------------------------------------

create-scan-table.R

 - *TO BE RUN LOCALLY*
 - Reads qc{qc_round}_review.xlsx created manually (see previous step) and generates next scan_table_for_registration.txt for failed scans only, using the suggested init_scale2 parameter. 
 - Export the table output back to the cluster, and rerun run-03.sh (example: rsync scan_table_for_registration2.txt rls83@login.hpc.cam.ac.uk:/home/rls83/rds/rds-rat-ex-vivo-QhAmE41a88w)
 - To rerun run-03.sh, you don't need to change anything in 03.align-to-atlas.R. Change the --array flag for number of scans to register (run cat scan_table_for_registration{qc_round}.txt | wc -l to check), change QC round in the -e and -o flag, and change qc_round in the script itself.
 - Repeat this process until all scans have been successfully registered (may take a few rounds of QC)


