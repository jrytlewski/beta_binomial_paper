## Differential Abundance Analysis of T cell clones from immunoSEQ(R) data
This repo contains data and python code associated with the Rytlewski et al manuscript on the beta-binomial model for differential abundance analysis of T cell clones.

### Healthy Donor Data
Immunosequencing data for the healthy donors is contained locally in this repo under healthy_donor_files.

### Urothelial Cancer Data
Immunosequencing data for the urothelial cancer patients can be freely downloaded from immuneACCESS: https://clients.adaptivebiotech.com/pub/snyder-2017-plosmedicine


## Usage
The following are example syntax that can be executed on the data found in the _healthy_donor_files_ folder. The configuration.ini settings should not be changed, unless specified, to reproduce results.


Example file for --batchFile:

```
subject_id  TSV 1  TSV 2
subjectA  fileA1.tsv.gz  fileA2.ts.gz
subjectB fileB1.tsv.gz  fileB2.tsv.gz
subjectC  fileC1.tsv.gz  fileC2.tsv.gz
```
Columns need to be labeled exactly as shown in above example. ".gz" should only be included when files are gzipped. The script is compatible with both .tsv and .tsv.gz immunoSEQ files.


#### Example 1 -- original binomial model
In the configuration.ini file, set _method = binomial_ before running.
```
python2.7 differential_abundance/rundiffabBatch_2017_09_24.py --batchfile differentialAbundance_batchfile.txt --config differential_abundance/configuration.ini --tsvDir healthy_donor_files/ --outDir diffab_results --parallel
```

#### Example 2 -- new beta binomial model
In the configuration.ini file, set _method = betabinomial_ before running.
```
python2.7 differential_abundance/rundiffabBatch_2017_09_24.py --batchfile differentialAbundance_batchfile.txt --config differential_abundance/configuration.ini --train differential_abundance/TrainingTSVs/replicates_Subject1_Standard.csv --tsvDir healthy_donor_files/ --outDir diffab_results --parallel
```

#### Disclaimer
For Research Use Only.

##### End of readme file.
