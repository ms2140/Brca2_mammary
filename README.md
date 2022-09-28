# Brca2_mammary
A transcriptional response to replication stress selectively expands a subset of BRCA2-mutant mammary epithelial cells

# Overview
This is the repository of scripts and data used for the paper "A transcriptional response to replication stress selectively expands a subset of BRCA2-mutant mammary epithelial cells"

# Documenation
Custom code was not generated for this manuscript. The following publicly available pipelines were utilized for analysis: raw sequencing data were processed using the 10X Genomics Cell Ranger pipeline (https://github.com/10XGenomics/cellranger), standard preprocessing and quality control of scRNA-seq data was based on Seurat Guided tutorials (http://satijalab.org/seurat/). 

# System Requirements
## Software version
R Version: 4.1.3 
Cell Ranger version: 6.0.2 
Seurat version: 4.0

## Hardware requirement 
Requires only a standard computer with enough RAM to support the in-memory operations

### R Dependencies
The image generation R-scripts depend on the following packages.

```
numpy
scipy
Cython
scikit-learn
pandas
seaborn
```

# Setting up the development environment:
- To build image and run from scratch:
  - Install [docker](https://docs.docker.com/install/)
  - Build the docker image, `docker build -t mgcpy:latest .`
    - This takes 10-15 mins to build
  - Launch the container to go into mgcpy's dev env, `docker run -it --rm --name mgcpy-env mgcpy:latest`
- Pull image from Dockerhub and run:
  - `docker pull tpsatish95/mgcpy:latest` or `docker pull tpsatish95/mgcpy:development`
  - `docker run -it --rm -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:latest` or `docker run -it --rm -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:development`


- To run demo notebooks (from within Docker):
  - `cd demos`
  - `jupyter notebook --ip 0.0.0.0 --no-browser --allow-root`
  - Then copy the url it generates, it looks something like this: `http://(0de284ecf0cd or 127.0.0.1):8888/?token=e5a2541812d85e20026b1d04983dc8380055f2d16c28a6ad`
  - Edit this: `(0de284ecf0cd or 127.0.0.1)` to: `127.0.0.1`, in the above link and open it in your browser
  - Then open `mgc.ipynb`

- To mount/load local files into docker container:
  - Do `docker run -it --rm -v <local_dir_path>:/root/workspace/ -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:latest`, replace `<local_dir_path>` with your local dir path.
  - Do `cd ../workspace` when you are inside the container to view the mounted files. The **mgcpy** package code will be in `/root/code` directory.


## MGC Algorithm's Flow
![MGCPY Flow](https://raw.githubusercontent.com/neurodata/mgcpy/master/MGCPY.png)

## Power Curves
- Recreated Figure 2 in https://arxiv.org/abs/1609.05148, with the addition of MDMR and Fast MGC
![Power Curves](https://raw.githubusercontent.com/neurodata/mgcpy/master/power_curves_dimensions.png)

# License

This project is covered under the **Apache 2.0 License**.
