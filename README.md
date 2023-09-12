
[Add PubMed link]: [![PubMed](https://salilab.org/imp-systems/static/images/pubmed.png)](https://pubmed.ncbi.nlm.nih.gov/36040254/)

[Add Zenodo link]: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6674232.svg)](https://doi.org/10.5281/zenodo.6674232)

# Integrative model of the WDR76-SPIN1-Nucleosome complex

This repository is of the integrative model of the WDR76-SPIN1-Nucleosome complex based on data from chemical crosslinking, X-ray crystallography, and structure prediction from [Alphafold](https://www.alphafold.ebi.ac.uk/entry/Q9H967). It contains input data, scripts for data preprocessing, modeling and results including bead models and localization probability density maps. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*).

These integrative structures will be deposited in the PDB-Dev database with accession codes ***AddPDBdev Ids***


***Keep or remove?***
![four_stage_figure](https://user-images.githubusercontent.com/8314735/165240737-c960c153-6af9-4014-9d76-32bc8c01cd03.png)


## Directory structure
1. [data](data/) : contains the subdirectories for the input data used for modeling all the subcomplexes.
2. [scripts](scripts/) : contains all the scripts used for modeling and analysis of the models.
3. [results](results/) : contains the models and the localization probability densities of the top cluster of the subcomplexes .
4. [test](test/) : scripts for testing the sampling


## Protocol
### Preprocessing the crosslinks
1. For crosslinks in sheetA of sheetA of `data/xlinks/original_suppmat_DataS3.xlsx`:  
    ```
    python get_protein_uniprot_mapping.py -x /home/shreyas/Dropbox/washburn_wdr_spin/xls_sheet1.csv
    ```
    Make a file `proteins_of_interest.txt` and use it to run:
    ```
    python get_protein_uniprot_mapping.py -x /home/shreyas/Dropbox/washburn_wdr_spin/xls_sheet1.csv -m mapping -p proteins_of_interest.txt
    ```
    Finally, to generate the input file for modeling, do:
    ```
    python xl_preprocessing.py ~/Dropbox/washburn_wdr_spin/xls_sheet1.csv uniprot_mapping.yaml
    ```
2. For crosslinks in sheetA of sheetA of `data/xlinks/original_suppmat_DataS3.xlsx`:  
    Run the following command to generate the crosslinks input file for modeling:
    ```
    python xl_change_protnames_from_ncbi2name.py
    ```


### Sampling
To run the sampling, run modeling scripts like this   
```
for runid in `seq 1 50` ; do mpirun -np 8 $IMP python scripts/modeling.py prod $runid ; done
```

where,   
`$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation).


### Analysis
#### 1. Getting the good scoring models
  Good scoring models were selected using `pmi_analysis` (Please refer to [pmi_analysis tutorial](https://github.com/salilab/PMI_analysis) for more detailed explaination) along with our `variable_filter_v1.py` script. These scripts are run as described below:
  1. First, run `run_analysis_trajectories_w_skip2.py` as follows:  
      `$IMP run_analysis_trajectories_w_skip2.py modeling run_ `  
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation),   
      `modeling` is the directory containing all the runs and   
      `run_` is the prefix for the names of individual run directories.  
      
  2. Then run `variable_filter_v1.py` on the major cluster obtained as follows:   
      `$IMP variable_filter_v1.py -c N -g MODEL_ANALYSIS_DIR`
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation),   
      `N` is the cluster number of the major cluster,   
      `MODEL_ANALYSIS_DIR` is the location of the directory containing the selected_models*.csv.   
      This can also be run using the `submit_variable_filter_v1.sh` script from the `scripts/analysis/pmi_analysis` directory.  
  _Please also refer to the comments in the `variable_filter_v1.py` for more details._

  3. The selected good scoring models were then extracted using `run_extract_models.py` as follows:   
      `$IMP python run_extract_models.py modeling run_ CLUSTER_NUMBER`   
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation),   
      `modeling` is the path to the directory containing all the individual runs and   
      `CLUSTER_NUMBER` is the number of the major cluster to be extracted.  
      
#### 2. Running the sampling exhaustiveness tests (Sampcon)
A separate directory named `sampcon` was created and a `density.txt` file was added to it. This file contains the details of the domains to be split for plotting the localisation probability densities. Finally, sampling exhaustiveness tests were performed using `imp-sampcon` as `$imp python $sampcon/pyext/src/exhaust.py -n wdr_spin -a -m cpu_omp -c 2 -d analysis/density.txt -gp -g 2.0  -sa ../model_analysis/A_gsm_clust1.txt -sb ../model_analysis/B_gsm_clust1.txt  -ra ../model_analysis/A_gsm_clust1.rmf3 -rb ../model_analysis/B_gsm_clust1.rmf3`. 

#### 3. Analysing the major cluster
1. Crosslink violations were analyzed as follows:   
    ```
    for xlfile in data/xlinks/modeling_xlfile_sheetA.dat data/xlinks/modeling_xlfile_sheetD.dat; do python get_xlviol_val_set_v2.py sampcon_cluster0_models.rmf3 $xlfile 35.0 & done
    ```   
      
2. Average distance maps were plotted for the models from the major cluster as follows:
    `scripts/analysis/cosmic_and_distance-maps/submit_contact_maps_all_pairs_surface.py`   
    This script calls the `scripts/analysis/cosmic_and_distance-maps/contact_maps_all_pairs_surface.py` script.
    _Please use `--help` for `contact_maps_all_pairs_surface.py` script for more details._

3. Plot distance versus model index plots for the models in the major cluster as follows:
    ```
    python scripts/analysis/plot_mdlike.py
    ```

### Results

For the simulations, the following files are in the [results](results/) directory
* `cluster_center_model.rmf3` : representative bead model of the major cluster
* `chimera_densities.py` : to view the localization densities (.mrc files)
* `xlviol` : Directory containing the logs for crosslink violations


## Information
**Author(s):** Xingyu Liu, Ying Zhang, Zihui Wen, Yan Hao, Charles Banks, Joseph Cesare, Saikat Bhattacharya, Shreyas Arvindekar, Jeffrey Lange, Brian Slaughter, Jay Unruh, Shruthi Viswanath, Laurence Florens, Jerry Workman, Michael Washburn  
**Date**: September 12th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** Not tested  
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:**  

