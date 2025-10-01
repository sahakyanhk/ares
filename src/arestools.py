import argparse
import os, re
import pandas as pd
import numpy as np
import shutil
import gzip
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors
import matplotlib.cm as cm

import RNA # pyright: ignore[reportMissingImports]

import warnings


parser = argparse.ArgumentParser(description="Analyse PFES")
parser.add_argument('-l', '--log', type=str, help='log file name', default='progress.log') 
parser.add_argument('-o', '--outdir', type=str, help='output directory name', default='visual_pfes_results')
parser.add_argument('-b', '--start', type=int, help='first point to read from trajectory', default=0)
parser.add_argument('-e', '--end', type=int, help='last point to read from trajectory', default=99999999)
parser.add_argument('--noplots', action='store_false', )


args = parser.parse_args()


#class VisualPFES():


def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(data, key=alphanum_key)

def extract_lineage(log):
    traj_len = len(log)
    #pop_size = len(log[log.gndx == 'gndx0'])

    print(f'Processing a trajectory with {traj_len} mutations')
    lineage = log.drop_duplicates('gndx').tail(1)
    df = lineage
    ndx = df.id.to_string(index=False)
        
    def return_ancestor(log, node):
        parent = log[log.id == node]
        parent = parent.drop_duplicates('sequence')
        return parent

    
    pbar = tqdm(desc='Extracting lineage')
    i=0
    while not df.empty:
        ndx = return_ancestor(log, ndx)
        df = ndx
        lineage = pd.concat([lineage, df], axis=0)
        ndx = ndx.prev_id.to_string(index=False)
        i=+1
        pbar.update(i)
    pbar.close()
    lineage = lineage.sort_index()
    ltail = lineage.tail(1)

    print(f"""
{ltail[['gndx',
        'id',
        'seq_len',
        'num_conts',
        'energy',
        'rfam_score',
        'rfam_id',
        'rfam_name',
        'score'
        ]].iloc[-1]}
{ltail.sequence.iloc[-1]}
{ltail.ss.iloc[-1]}
""")
    return lineage

def extract_sequences(log):
    fasta  = "fasta"
    return fasta

#======================= make separate plots =======================#
def make_plots(log, bestlog, lineage):

    ms=0.1
    lw=1.4
    dpi=500

    os.makedirs(plotdir, exist_ok=True)
    for colname in log.keys(): 
        if not colname in ['gndx', 'id', 'seq_len', 'num_conts', 'sel_mode', 'rfam_id', 'rfam_name', 'sequence', 'mutation', 'prev_id', 'ss']:
                fig, ax1 = plt.subplots(figsize=(9, 3))
                ax1.plot(log[colname],'.', markersize=ms,    color='silver', label='all mutations')
                ax1.plot(bestlog[colname],'-', linewidth=lw, label='best of the generation')
                ax1.plot(lineage[colname],'-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={len(lineage[colname])})')
                ax1.legend(loc ="lower right")
                ax1.grid(True, which="both",linestyle='--', linewidth=0.3)
                ax1.set(xlabel="Total number of mutations", ylabel=colname.capitalize())
                #ax2 = ax1.twiny()
                #ax2.plot(lineage[colname].tolist(),'-', linewidth=lw, color='mediumslateblue')
                #ax2.set(xlabel="Lineage length")
                fig.tight_layout()
                fig.savefig(plotdir + colname + '.png', dpi=dpi)
                fig.clf()

#======================= Summary plot =======================#
def make_summary_plot(log, bestlog, lineage):
    
    ms=0.1
    lw=1.0
    dpi=500

    fig, axs = plt.subplots(3,2, figsize=(10, 8))

    fig.suptitle("")

    L = len(lineage)
    axs[0,0].plot(log.energy, '.', markersize=ms,    color='silver', label='all mutations')
    axs[0,0].plot(bestlog.energy, '-', linewidth=lw, label='best of the generation')
    axs[0,0].plot(lineage.energy, '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[0,0].set(xlabel=None, ylabel='SS energy')
    axs[0,0].grid(True, which="both",linestyle='--', linewidth=0.5)
    axs[0,0].set_xticklabels([])

    axs[1,0].plot(log.rfam_score, '.', markersize=ms,    color='silver', label='all mutations')
    axs[1,0].plot(bestlog.rfam_score, '-', linewidth=lw, label='best of the generation')
    axs[1,0].plot(lineage.rfam_score, '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[1,0].set(xlabel=None, ylabel='RFAM score')
    axs[1,0].grid(True, which="both",linestyle='--', linewidth=0.5)
    axs[1,0].set_xticklabels([])

    axs[2,0].plot(log.score,  '.', markersize=ms,    color='silver', label='all mutations')
    axs[2,0].plot(bestlog.score, '-', linewidth=lw,  label='best of the generation')
    axs[2,0].plot(lineage.score,  '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[2,0].set(xlabel='Total number of mutations', ylabel='Score')
    axs[2,0].grid(True, which="both",linestyle='--', linewidth=0.5)
    axs[2,0].legend(loc ="lower right")

    axs[0,1].plot(log.seq_len, '.', markersize=ms,    color='silver', label='all mutations')
    axs[0,1].plot(bestlog.seq_len, '-', linewidth=lw, label='best of the generation')
    axs[0,1].plot(lineage.seq_len, '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[0,1].set(xlabel=None, ylabel='Sequence length')
    axs[0,1].grid(True, which="both",linestyle='--', linewidth=0.5)
    axs[0,1].set_xticklabels([])

    axs[1,1].plot(log.gc_cont, '.', markersize=ms,    color='silver', label='all mutations')
    axs[1,1].plot(bestlog.gc_cont, '-', linewidth=lw, label='best of the generation')
    axs[1,1].plot(lineage.gc_cont, '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[1,1].set(xlabel=None, ylabel='GC content')
    axs[1,1].grid(True, which="both",linestyle='--', linewidth=0.5)
    axs[1,1].set_xticklabels([])

    axs[2,1].plot(log.num_conts, '.', markersize=ms,  color='silver', label='all mutations')
    axs[2,1].plot(bestlog.num_conts, '-', linewidth=lw, label='best of the generation')
    axs[2,1].plot(lineage.num_conts, '-', linewidth=lw, color='mediumslateblue', label=f'lineage (L={L})')
    axs[2,1].set(xlabel='Total number of mutations', ylabel='Number of contacts')
    axs[2,1].grid(True, which="both",linestyle='--', linewidth=0.5)

    #plt.xticks(rotation=45)

    #for ax in axs.flat:
    #   ax.set(xlabel='x-label', ylabel='y-label')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    #for ax in axs.flat:
    #   ax.label_outer()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir,'Summary.png'), dpi=dpi)

#======================= seconday structure plot =======================#

def make_ss_plot(lineage):

    dpi=500

    max_seq_len = int(max(lineage.seq_len))
    lineage_len = len(lineage)
    ss = list(lineage.ss)

    def parse_stems(dotbracket):
        """Return stems as lists of paired indices [(i,j), ...]"""
        stack = []
        pairs = []
        for i, c in enumerate(dotbracket):
            if c == "(":
                stack.append(i)
            elif c == ")":
                j = stack.pop()
                pairs.append((j, i))
        pairs.sort()

        # group into stems (consecutive nested pairs)
        stems = []
        current = [pairs[0]]
        for (i,j), (i2,j2) in zip(pairs, pairs[1:]):
            if i2 == i+1 and j2 == j-1:  # stacked
                current.append((i2,j2))
            else:
                stems.append(current)
                current = [(i2,j2)]
        stems.append(current)
        return stems

    fig, ax = plt.subplots(figsize=(10, len(ss)*0.1))


    for y, struct in enumerate(ss):
        stems = parse_stems(struct)
        colors = cm.tab20(np.linspace(0,1,len(stems)))  # type: ignore # distinct colors
        color_map = {}

        # assign color to each position in stems
        for c, stem in zip(colors, stems):
            for (i,j) in stem:
                color_map[i] = c
                color_map[j] = c

        # plot structure
        for x, c in enumerate(struct):
            if c == ".":
                ax.text(x, y, ".", color="gray", fontsize=7, ha="center", va="center")
            elif c == "(":
                ax.text(x, y, "<", color=color_map.get(x,"blue"), fontsize=7, ha="center", va="center")
            elif c == ")":
                ax.text(x, y, ">", color=color_map.get(x,"red"), fontsize=7, ha="center", va="center")


    ax.set_xlim(-1, len(ss[0]))
    ax.set_ylim(-1, len(ss))
    ax.invert_yaxis()
    ax.set_xlabel("Nucleotide position")
    ax.set_ylabel("Structure index")
    ax.set_title("RNA secondary structures (colored stems)")
    ax.axis("tight")
    plt.savefig(os.path.join(outdir,'Secondary_structures.png'), dpi=dpi) 


#======================= seconday structure movie =======================#
import cairosvg
from moviepy import ImageClip, concatenate_videoclips

def meke_ss_movie(lineage): 
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    
    os.makedirs(tmp_svg, exist_ok=True)
    os.makedirs(tmp_png, exist_ok=True)

    sequences = list(lineage.sequence)
    secondary_structres = list(lineage.ss)

    assert len(sequences) == len(secondary_structres)

    i=0
    for seq, ss in zip(sequences, secondary_structres): 
        
        RNA.plot_structure_svg(f"{tmp_svg}rna{i}.svg", seq, ss)
        cairosvg.svg2png(url=f"{tmp_svg}rna{i}.svg", write_to=f"{tmp_png}rna{i}.png")
    
        i+=1

    
    clips = [ImageClip(tmp_png + m).with_duration(0.1)
           for m in sorted_alphanumeric(os.listdir(tmp_png))]
             
    concat_clip = concatenate_videoclips(clips, 
                                        method="compose", 
                                        bg_color=(255, 255, 255))
    concat_clip.write_videofile(f'{outdir}/rne_evolution_movie.mp4', 24)


#======================= functions end here =======================#

outdir = args.outdir 
os.makedirs(outdir, exist_ok=True)
plotdir = os.path.join(outdir, 'plots/')
tmp = os.path.join(outdir, 'tmp/')
tmp_svg = os.path.join(outdir, 'tmp/svg/')
tmp_png = os.path.join(outdir, 'tmp/png/')

log = pd.read_csv(args.log, sep='\t', comment='#')
log = log.iloc[args.start:args.end]

bestlog = log.groupby('gndx').head(1)
bestlog.to_csv(os.path.join(outdir, 'bestlog.tsv'), sep='\t', index=False, header=True)


print('==================================')
lineage = extract_lineage(log)
lineage.to_csv(os.path.join(outdir, 'lineage.tsv'), sep='\t', index=False, header=True)





if args.noplots:

    # print('making summary plot')
    # make_summary_plot(log, bestlog, lineage)

    # print('making plots')
    # make_plots(log, bestlog, lineage)

    # print('making secondary structure plot')
    # make_ss_plot(lineage)
    
    print('making secondary structure movie')
    meke_ss_movie(lineage)

print('=================================='"\n")

