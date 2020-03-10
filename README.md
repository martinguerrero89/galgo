#  Galgo: An evolutionary machine learning algorithm for the identification and study of prognostic gene expression signatures in cancer

**Motivation**: Clustering analysis has been long used to find underlying structures in different omics data such as gene expression profiles. This data typically presents high number of dimensions and has been used successfully to find co-expressed genes in samples that share similar molecular and clinical characteristics. Nevertheless, the clustering results are highly dependent of the features used and the number of clusters considered, while the partition obtained does not guarantee clinically relevant findings.

**Methods**: We propose a multi-objective optimization algorithm for disease subtype discovery based on a non-dominated sorting genetic algorithm. Our proposed framework combines the advantages of clustering algorithms for grouping heterogeneous omics data and the searching properties of genetic algorithms for feature selection and optimal number of clusters determination to find features that maximize the survival difference between subtypes while keeping cluster consistency high.

### To run the algorithm all the data must be in the *Data* folder

* Open R and set the main Galgo folder as working directory
* First, be sure you have all the needed libraries (listed in the 'libraries.R' file) installed and properly configured.
You can install the released version of galgo using devtools with:

``` r
devtools::install_github("https://github.com/harpomaxx/galgo")
```

* Of importance is to have `gpuR` package propperly installed and working. To instruction on the installation process and troubleshoot follow: https://github.com/cdeterman/gpuR/wiki
By default galgo runs some portions of its code in GPU, provided by the gpuR package. Before installing gpuR, the opencl backend should be configured. 

In linux systems install lastest nvidia cuda drivers and the opencl backend.

```
       apt-get install nvidia-418 nvidia-opencl-icd-418 libcuda1-418
       apt-get install opencl-headers  ocl-icd-opencl-dev
       
```

* Once the repository is cloned and all the dependencies are install run:

```
      source("./libraries.R")
      
```

* For each cancer type run the following set of commands

##Breast cancer
```
      #Create data
      source("./Data/BRCA_data.R")
      
      #Run galgoR
      source("./Galgo_BRCA.R")
      
```
##Colorectal cancer
```
      #Create data
      source("./Data/CRC_data.R")
      
      #Run galgoR
      source("./Galgo_CRC.R")
      
```
##Lung adenocarcinoma
```
      #Create data
      source("./Data/LUAD_data.R")
      
      #Run galgoR
      source("./Galgo_LUAD.R")
      
```

##High grade serous ovarian cancer
```
      #Create data
      source("./Data/HGSOC_data.R")
      
      #Run galgoR
      source("./Galgo_HGSOC.R")
      
```
