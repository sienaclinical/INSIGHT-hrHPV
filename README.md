## Integrative Analyses of Individual Patient Genomic-Data to Discover Novel Biomarkers

Project code name: INSIGHT

### About the project
Despite advancements in its affordability, sequencing analyses require substantial financial and computational resources, limiting widespread implementation, especially in resource-limited settings. Individual studies often discover promising biomarkers, yet they reflect local optima with limited generalization. Integrative analysis of multiple datasets across diverse populations has potential to improve understanding of disease progression and treatment outcomes. The project utilizes data aggregation techniques to use publicly available individual patients genomic-data sets for cancer prediction.

### Use case on Cervical pre-cancer on hrHPV women population
Uncovering global optimum biomarkers capable of distinguishing between low-grade squamous intraepithelial lesions (LSIL, refer to as “control” group) and high-grade squamous intraepithelial lesions (HSIL, refer to as “case” group) in women who test positive for high-risk HPV. The project focuses on integrating DNA sequencing datasets from previously published studies, as DNA sequencing is widely utilized in the discovery phase for identifying potential biomarkers.

### Datasets
Search strategies on GEO repository (not elaborated here) yields in four previously published studies that met requirements for the aforementioned use case.
1. GSE99511 - Raw datasets from Illumina HumanMethylation450 BeadChip platform
2. GSE143752 - Preprocessed data from Infinium MethylationEPIC platform
3. GSE186835 - Raw datasets from Illumina HumanMethylation450 BeadChip platform 
4. GSE287994 - Raw datasets from Infinium MethylationEPIC platform

Total number of samples
- Control: 171 samples
- Case: 216 samples

### Workflow
#### 1️⃣ Preprocessing datasets
script: 

#### 2️⃣ Calculate "effect size" on individual study
#### 3️⃣ Aggregating multiple studies via MA-approach
