#!/bin/bash
sbatch <<EOF
#!/bin/sh
#SBATCH --job-name=$(basename $1)
#SBATCH --cpus-per-task=12
#SBATCH --mem=50g
#SBATCH --partition=norm
#SBATCH --time=12:00:00
#SBATCH --output="${1}.out"


python src/ares.py --nobackup \
                    -iseq random \
                    --random_seq_len 90 \
                    -ps 100 -ng 500 -ed rna -em single_chain \
                    -b 1 -o "$1"

python src/arestools.py -l "$1"/progress.log -o "$1"


EOF
