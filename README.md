## A network medicine approach to investigation and population-based validation of disease manifestations and drug repurposing for COVID-19

### Network proximity code
The proximity code can be found in the `proximity` folder. First, you need to decompress the file `HumanInteractome.7z` in the same folder.

* HumanInteractome.tsv - the human interactome, including sources and evidence types
* HumanInteractome.npy - numpy matrix of all precalculated shortest distance in the interactome
* DrugTargetNetwork.txt - drug target network
* network_proximity.py - code to compute closest network proximity. See below

The `network_proximity.py` supports two modes. To run the program (Python 3), numpy and networkx need to be installed.
```
pip install numpy networkx
```

#### To compute the closest proximity between two gene lists
```
python network_proximity.py path_to_gene_list_1 path_to_gene_list_2 number_of_repeat random_seed
```
Example
```
python network_proximity.py example/Asthma.txt example/SARS2-DEG_lung.txt 1000 11096
```
#### To screen all 3000 drugs for a single gene list
```
python network_proximity.py DRUG path_to_gene_list number_of_repeat random_seed
```
Example
```
python network_proximity.py DRUG example/SARS2-DEG.txt 1000 1024
```

Please see https://github.com/ChengF-Lab/GPSnet and https://chengf-lab.github.io/PDGPS/ for more details and explanations

### Supplemental Files
Supplemental figures and tables in found in supplemental_files/

* S1 Table. Summary of the data sets used in this study.
* S2 Table. Five SARS-CoV-2 target data sets used in this study.
* S3 Table. Additional virus target lists for comparisons with SARS-CoV-2 targets.
* S4 Table. Disease-associated genes.
* S5 Table. COVID-19 clinical studies used in the meta-analysis.
* S6 Table. Network proximity results for 2,938 drugs against the SARS-CoV-2 data sets.
* S7 Table. Proposed repurposable drugs and their antiviral profiles.
* S1 File. Network file for the global network of disease manifestations associated with human coronavirus.
