
python ares.py --nobackup -iseq randoms --random_seq_len 50  -o test -ps 100 -ng 100 -ed rna -em single_chain -b 0.5

python arestools.py -l test/progress.log -o test


RNAplot