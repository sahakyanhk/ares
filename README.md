# ARES: aproximated RNA evolution simulator


### Installation

Install dependencies:

[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA?tab=readme-ov-file#installation)\
`conda install viennarna`

[Infernal](http://eddylab.org/infernal)\
`conda install -c bioconda infernal`

Download code and prepare Rfam database
```
wget https://github.com/sahakyanhk/ares/archive/refs/heads/main.zip -O ares.zip; unzip ares.zip

bash ares/src/make_rfam_db.sh

python src/ares.py
```

### Usage

```
python src/ares.py --nobackup -iseq random --random_seq_len 90  -o test/run1 -ps 100 -ng 500 -ed rna -em single_chain -b 1

python src/arestools.py -l test/run1/progress.log -o test/run1
```
or submit a sbatch job with SLURM:

```
bash run_ares.sbatch test/run
```