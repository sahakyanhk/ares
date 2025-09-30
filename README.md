# ARES: aproximated RNA evolution simulator


### Installation

Install ependencies

(ViennaRNA)[https://github.com/ViennaRNA/ViennaRNA?tab=readme-ov-file#installation]
conda install viennarna

(Infernal)[[http://eddylab.org/infernal]]
conda install -c bioconda infernal

Downlaod code and prepare Rfam database
```
wget https://github.com/sahakyanhk/ares/archive/refs/heads/main.zip -O ares.zip; unzip ares.zip

bash ares/src/make_rfam_db

python src/ares.py
```
