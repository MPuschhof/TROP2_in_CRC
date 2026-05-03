# TROP2 in CRC

**TROP2 targeting reveals therapy-driven cell state transitions in human colorectal cancer **

Metastatic colorectal cancer (CRC) remains a leading cause of cancer-related mortality and is driven by pronounced tumour cell plasticity. Here we identify the transmembrane glycoprotein trophoblast cell-surface antigen 2 (TROP2) as a marker of poor-prognosis CRC associated with WNT-low, fetal-like tumour cell states linked to metastasis and therapy resistance. Functional analyses demonstrate that TROP2⁺ cells exhibit context-dependent stem-like capacity and the ability to initiate metastatic outgrowth. Given that these detrimental tumour states converge on the cell-surface antigen TROP2, we explored therapeutic targeting of this cell population using clinically relevant TROP2-directed antibody-drug conjugates (ADCs). Time-resolved analyses reveal therapy-associated dynamics in tumour cell state composition between WNT-high LGR5⁺ intestinal stem cell programs and WNT-low TROP2⁺ fetal-like states. Conventional chemotherapy promotes the induction of TROP2-expressing cells, whereas TROP2-ADCs selectively target these populations and remodel the tumour cell-state landscape. Exploiting this plasticity, combined chemotherapy and TROP2 targeting enhances antitumour efficacy in patient-derived models. Together, our findings identify TROP2 as a therapeutic vulnerability of poor-prognosis CRC and highlight the importance of targeting tumour cell states to improve therapeutic efficacy and overcome adaptive resistance in advanced disease.

## Data access

Raw and processed data objects have been made available here through arrayExpress / European Nucleotide Archive (ENA) under the following accession numbers: CRC cohort tissue-seq [E-MTAB-16583](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16583?query=E-MTAB-16583) and PDOX-seq [E-MTAB-16585](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16585?query=E-MTAB-16585); SG treatment of subcutaneous tumours [E-MTAB-16433](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16433?query=E-MTAB-16433); time course of SG treatment of liver metastasis [E-MTAB-16849](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16849?query=E-MTAB-16849), time course of SG treatment of PDOs in vitro [E-MTAB-16843](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16843?query=E-MTAB-16843); time course of FOLFIRI treatment of PDOs [E-MTAB-16836](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16836?query=E-MTAB-16836); MDO-derived tumors VAKPS and VKPN [E-MTAB-16835](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-16835?query=E-MTAB-16835).

The data can be downloaded and re-assembled to anndata files as demonstrated in the notebook [get_data_from_arrayExpress.ipynb](https://github.com/ORCA-HD/TROP2_in_CRC/blob/main/figures/get_data_from_arrayExpress.ipynb).

### Expected runtime

Depending on the size of the dataset and the download speed, the time required for downloading and re-assembling a dataset can be quite variable. For us, it’s typically done within 10 minutes.


## Running the analysis

Having downloaded and re-assembled the processed data as indicated in the section on [data access](https://github.com/ORCA-HD/TROP2_in_CRC/tree/main#data-access), the analysis can be reproduced using the notebooks in [figures](https://github.com/ORCA-HD/TROP2_in_CRC/tree/main/figures). The analysis has been performed using different micromamba environments. 

The environments can be recreated using the specifics listed in [env_specs](https://github.com/ORCA-HD/TROP2_in_CRC/tree/main/env_specs). For instance, use `micromamba create --file env_general.yml` to create the general analysis environment.

Most analyses were performed using the general environment *env_general*. Analysis requiring trajectory inference through palantir were done with *env_palantir* and pseudobulk analysis with the previous version of decoupler (1.8.0) were run with *env_decoupler_v1.8*.

In addition, the session information is provided at the end of each notebook.

### Expected runtime

Recreating the analysis from the processed data is done within minutes for each notebook as the results of memory-heavy computation steps has been deposited to arrayExpress. Where required, loading magic-imputed counts can take several minutes on its own and might require up to 32 GB of memory. 


## Publication

Vaquero, Georgakopoulos, Puschhof *et al.* (submitted, 2026). TROP2 targeting reveals therapy-driven cell state dynamics in colorectal cancer.
