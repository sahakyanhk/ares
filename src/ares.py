import os
import sys
import shutil
import pandas as pd
import argparse
import typing as T
import threading
import gzip
import time
from datetime import datetime

from evolution import Evolver
from rnautils import rnafold, rna_seq_search

def backup_output(outpath):
    print(f'\nSaving output files to {args.outpath}')
    if os.path.isdir(outpath): 
        backup_list = []
        last_backup = int()
        for dir_name in os.listdir():
            if dir_name.startswith(outpath + '.'):
                backup=(dir_name.split('.')[-1])
                if backup.isdigit(): 
                    backup_list.append(backup)
                    last_backup = int(max(backup_list))
        print(f'\n{outpath} already exists, renameing it to {outpath}.{str(last_backup +  1)}') 
        os.replace(outpath, outpath + '.' + str(last_backup +  1))

def create_batched_rna_sequence_dataset(sequences: T.List[T.Tuple[str, str]], max_seq_per_batch: int = 50
) -> T.Generator[T.Tuple[T.List[str], T.List[str]], None, None]:
    batch_headers, batch_sequences, num_sequences= [], [], 0 
    for header, seq in sequences:
        if num_sequences > max_seq_per_batch:
            yield batch_headers, batch_sequences
            batch_headers, batch_sequences, num_sequences= [], [], 0

        batch_headers.append(header)
        batch_sequences.append(seq)
        num_sequences += 1
        if num_sequences > args.pop_size: #TODO test this with args.pop_size / 4 and lartge pop size
           yield batch_headers, batch_sequences
           batch_headers, batch_sequences, num_sequences= [], [], 0
    yield batch_headers, batch_sequences


def extract_results(gen_i, headers, sequences, secondary_structures, energies, rfam_score_list, rfam_id_list, target_name_list) -> None:
    global new_gen #this will be modified in the fold_evolver()

    for meta_id, seq, ss, energy, rfam_score, rfam_id, rfam_name in \
        zip(headers, sequences, secondary_structures, energies, rfam_score_list, rfam_id_list, target_name_list):
        
        all_seqs = seq.split(':')
        seq = all_seqs[0]
        seq_len = len(seq)
        
        id_data = meta_id.split('_')

        id = id_data[0]
        prev_id = id_data[1]
        mutation = id_data[2]

        #=======================================================================# 
        #================================SCORING================================# 

        gc_cont = round(((seq.count('C') + seq.count('G')) / seq_len )* 100, 2)
        num_conts = ss.count("(")
        #score = num_conts / seq_len * energy
        score = (rfam_score - 0.5 * energy) / seq_len
        #================================SCORING================================#
        #=======================================================================# 

        iterlog = pd.DataFrame({'gndx': gen_i,
                                'id': id, 
                                'seq_len': seq_len,
                                'num_conts': num_conts, 
                                'sel_mode': args.selection_mode,
                                'beta': args.beta,
                                'energy': round(energy, 3), 
                                'rfam_score': rfam_score,
                                'rfam_id': rfam_id,
                                'rfam_name': rfam_name,
                                'gc_cont': gc_cont,
                                'score': round(score, 3),
                                'sequence': seq, 
                                'mutation': mutation,
                                'prev_id': prev_id,
                                'ss': ss,
                                }, index=[0])
        
        new_gen = pd.concat([new_gen, iterlog], axis=0, ignore_index=True) 

    
def fold_evolver(args, evolver, logheader, init_gen) -> None: 

    tmp = args.outpath + "/tmp"

    os.makedirs(args.outpath, exist_ok=True)
    with open(os.path.join(args.outpath, args.log), 'w') as f:
        f.write(logheader)

    
    #creare an initial pool of sequences with pop_size
    columns=['gndx',
              'id',
              'seq_len',
              'num_conts',
              'sel_mode',
              'beta',
              'energy',
              'rfam_score',
              'rfam_id',
              'rfam_name',
              'gc_cont',
              'score',
              'sequence',
              'mutation',
              'prev_id',
              'ss']

    ancestral_memory = pd.DataFrame(columns=columns)
    ancestral_memory.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=True, sep='\t') #write header of the progress log
    
    #mutate seqs from init_gen and select the best N seqs for the next generation    
    for gen_i in range(args.num_generations):
        n = 0
        global new_gen #this will be modified in the extract_results() 
        new_gen = pd.DataFrame(columns=columns)
        #now = datetime.now()
        generated_sequences = []
        mutation_collection = []

        for prev_id, sequence in zip(init_gen.id, init_gen.sequence):
            seq, mutation_data = evolver.mutate(sequence)
            
            #check if the mutated seqeuece was already predicted
            seqmask = ancestral_memory.sequence == seq 
            
            #if --norepeat and seq is in the ancestral_memory mutate it again
            if args.norepeat and seqmask.any():  
                while seqmask.any():
                    seq, mutation_data = evolver.mutate(seq)
                    seqmask = ancestral_memory.sequence == seq 

            id = f"g{gen_i}seq{n}_{prev_id}_{mutation_data}"; n+=1 # gives an unique id even if the same sequence already exists            

            if seqmask.any(): #if sequence already exits do not predict a structure again 
                repeat = ancestral_memory[seqmask].drop_duplicates(subset=['sequence'], keep='last') 
                new_gen = pd.concat([new_gen, repeat])
            else:
                generated_sequences.append((id, seq)) 
                mutation_collection.append(mutation_data) 
        
        batched_sequences = create_batched_rna_sequence_dataset(generated_sequences, args.max_seq_per_batch)

        #predict data for the new batch
        for headers, sequences in batched_sequences:
            secondary_structures, energies  = rnafold(sequences) # type: ignore
            rfam_score_list, rfam_id_list, target_name_list = rna_seq_search(headers, sequences, tmp) # data is sorted in the same order as headers

            #run extract_results() in beckground and imediately start next the round of model.infer()
            trd = threading.Thread(target=extract_results, \
                        args=(gen_i, headers, sequences, secondary_structures, energies, rfam_score_list, rfam_id_list, target_name_list))
            trd.start()

        print(new_gen.tail(args.pop_size).drop(columns=['gndx','ss'],
                                            axis=1).to_string(index=False, header=False))


        #print(f"#GENtime {datetime.now() - now}")
        ancestral_memory =  pd.concat([ancestral_memory, init_gen])

        #select the next generation 
        init_gen = evolver.select(new_gen, init_gen, args.pop_size, args.selection_mode, args.norepeat, args.beta)
        init_gen.gndx = f'gndx{gen_i}' #assign a new gen index
        init_gen.to_csv(os.path.join(args.outpath, args.log), mode='a', index=False, header=False, sep='\t')

 
#================================FOLD_EVOLVER================================#
#============================================================================#


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Sample sequences based on a given structure.'
    )
    parser.add_argument(
            '-em', '--evolution_mode', type=str,
            help='evolution mode: single_chain, inter_chain, multimer',
            default='single_chain',
    )
    parser.add_argument(
            '-sm', '--selection_mode', type=str,
            help='selection mode\n options: strong, weak ',
            default="weak"
    )
    parser.add_argument(
            '-b', '--beta', type=float,
            help='selection strength',
            default=1,
    )
    parser.add_argument(
            '-iseq', '--initial_seq', type=str,
            help='a sequence to initiate with, if "random" pop_size random sequences will ' \
            'be generated, the length of the random sequences can be assigned with "--random_seq_len"',
            default='random'
    )
    parser.add_argument(
            '-l', '--log', type=str,
            help='log output',
            default='progress.log',
    )   
    parser.add_argument(
            '-o' ,'--outpath', type=str,
            help='output filepath for saving sampled sequences',
            default='output',
    )
    parser.add_argument(
            '-ng', '--num_generations', type=int,
            help='number of generations',
            default=100,
    )
    parser.add_argument(
            '-ps', '--pop_size', type=int,
            help='population size',
            default=10,
    )
    parser.add_argument(
            '--max_seq_per_batch', type=int,
            help='population size',
            default=50,
    )
    parser.add_argument(
            '-ed', '--evoldict', type=str,
            help='population size',
            default='flatrates',
    )
    parser.add_argument(
            '--random_seq_len', type=int,
            help='a sequence to initiate with',
            default=24,
    )
    parser.add_argument(                      
            '--norepeat', action='store_true', 
            help='do not generate and/or select the same sequences more than once', 
    )
    parser.add_argument(
            '--nobackup', action='store_true', 
            help='overwrite files if exists',
    )


    args = parser.parse_args()
    evolver = Evolver(args.evoldict)

    now = datetime.now() # current date and time
    date_now = now.strftime("%d-%b-%Y")
    time_now = now.strftime("%H:%M:%S")
    


    logheader = f'''#======================== PFESv0.1 ========================#
#====================== {date_now} =======================#
#======================== {time_now} ========================#
#WD: {os.getcwd()}
#$pfes.py {' '.join(sys.argv[1:])}
#
#====================  pfes input params ==================#
#
#--evolution_mode, -em \t\t = {args.evolution_mode}
#--selection_mode, -sm\t\t = {args.selection_mode}
#--initial_seq, -iseq\t\t = {args.initial_seq}
#--pop_size, -ps\t\t = {args.pop_size}
#--evoldict, -ed\t\t = {args.evoldict}
#--log, -l\t\t\t = {args.log}
#--outpath, -o\t\t\t = {args.outpath}
#--random_seq_len\t\t = {args.random_seq_len}
#--beta, -b\t\t\t = {args.beta}
#--norepeat\t\t\t = {args.norepeat}
#--nobackup\t\t\t = {args.nobackup}
# evolution dictionary = {evolver.evoldict}
# evolution dictionary normalized = {evolver.evoldict_normal}
#==========================================================#
'''
    
    print(logheader)

    #backup if output directory exists
    if args.nobackup:
        if os.path.isdir(args.outpath):
            print(f'\nWARNING! Directory {args.outpath} exists, it will be replaced!')
            shutil.rmtree(args.outpath)
        os.makedirs(args.outpath)
    else:
        backup_output(args.outpath)


    #create the initial generation
    if args.initial_seq == 'random':
        randomsequence = evolver.randomseq(args.random_seq_len)
        init_gen = pd.DataFrame({'id': ['initseq'] * args.pop_size, 
                                 'sequence': [randomsequence] * args.pop_size,
                                 'score': [1e-99] * args.pop_size})
    elif args.initial_seq == 'randoms':
        init_gen = pd.DataFrame({'id': [f'initseq{i}' for i in range(args.pop_size)], 
                                 'sequence': [evolver.randomseq(args.random_seq_len) for i in range(args.pop_size)],
                                 'score': [1e-99] * args.pop_size})
    #elif args.initial_seq == 'c':
    #    init_gen = pd.read_csv('test.chk', sep='\t')
    else: 
        init_gen = pd.DataFrame({'id': ['initseq'] * args.pop_size, 
                                 'sequence': [args.initial_seq] * args.pop_size,
                                 'score': [1e-99] * args.pop_size})
        
    init_gen["ss"], init_gen["energy"] = rnafold(list(init_gen.sequence))   # type: ignore
    init_gen["seq_len"] = [len(seq) for seq in init_gen.sequence]
    init_gen["num_conts"] = [ss.count("(") for ss in init_gen.ss]


    print('running RFES... \n')
    if args.evolution_mode == "single_chain":
        fold_evolver(args, evolver, logheader, init_gen)
    elif not args.evolution_mode in ['single_chain', 'inter_chain', 'multimer']:
        print("Unknown PFES mode: aveilable options are: single_chain, inter_chain or multimer")





