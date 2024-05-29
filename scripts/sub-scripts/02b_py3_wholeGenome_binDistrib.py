import sys

########################  LOAD INPUT VARIABLES  ##########################
# Path to Sample table
sample_table = sys.argv[1]
# File with gathered results from alignment and peak calling
gatheredResults = sys.argv[2]
# Path with "sps_id" folders containing bam files from each species
inPath = sys.argv[3]
# Path to output PDFs
out_dir = sys.argv[4]
# Maximum number of CPUs to use
nCPU = sys.argv[5]
# Path to iChIPv2 Conda env bin folder
condaEnv = sys.argv[6]

# Other variables
# Resolution (in bp) of the computed bigwig files
resolution = 1000
# Set to True if you want to store output plots
saveFig = True



#######################################  Libraries  ######################################
import pyBigWig
import os
import subprocess
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
from collections import defaultdict
import re
import pandas as pd
import pysam
from scipy.stats import gaussian_kde
from scipy.stats import kde as kde


def runCommand(cmd):
    # run job
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    

def gaussAreaRight(x0, x1, ax):
    col1 = sns.color_palette()[0]
    col2 = sns.color_palette()[1]

    #x0 = df_start.loc[pos, celchi]
    #x1 = df_start.loc[pos == False, celchi]

    kde0 = gaussian_kde(x0, bw_method=0.3)
    kde1 = gaussian_kde(x1, bw_method=0.3)

    xmin = min(x0.min(), x1.min())
    xmax = max(x0.max(), x1.max())
    dx = 0.2 * (xmax - xmin) # add a 20% margin, as the kde is wider than the data
    xmin -= dx
    xmax += dx

    x = np.linspace(xmin, xmax, 500)
    kde0_x = kde0(x)
    kde1_x = kde1(x)
    inters_x = np.minimum(kde0_x, kde1_x)

    ax.plot(x, kde0_x, color=col1, label='Control')
    ax.fill_between(x, kde0_x, 0, color=col1, alpha=0.2)
    ax.plot(x, kde1_x, color=col2, label='ChIP')
    ax.fill_between(x, kde1_x, 0, color=col2, alpha=0.2)

    area_inters_x = np.trapz(inters_x, x)

    medX = np.median(x0)

    # Add vertical line in cutting points
    idx = np.argwhere(np.diff(np.sign(kde0_x - kde1_x))).flatten()
    cutP1 = None
    cutP2 = None
    #print(idx)
    for cut in idx:
        if x[cut] > medX and cutP1 == None:
            #print(x[cut])
            cutP1 = cut
        elif (x[cut] > medX) and (cutP2 == None):
            cutP2 = cut
    if cutP2 == None:
        cutP2 = len(x)-1
    for cut in [cutP1, cutP2]:
        ax.axvline(x[cut], linestyle='--', color='grey')

    # area under the lower kde, from the second intersection point to the last rightmost point
    inters_2 = kde1_x[cutP1:] - kde0_x[cutP1:]
    area2 = np.trapz(inters_2, x[cutP1:]) 
    ax.plot(x[cutP1:], inters_2, color='r')
    ax.fill_between(x[cutP1:], inters_2, 0, facecolor='none', 
                     edgecolor='r', hatch='xx', label='AreaPass')

    handles, labels = ax.get_legend_handles_labels()
    labels[2] += f': {area2 * 100:.1f} %'
    ax.legend(handles, labels, title='Location')
    ax.set_title('Signal vs IgG log2(FC + 1)')
    #ax.tight_layout()

    return round(area2 * 100, 2)
#--------------------------------------------------------------
## VARIABLES NOT TO BE CHANGED
# Output files path
outdata = f"{out_dir}/temp"
# Output plot files main path
outplot = f"{out_dir}/QC"

# Bigwig output directory
outBigW= f"{outdata}/allBigwig_{resolution//1000}kb"

# Color ranges fro FRiP scores
fripColors = {(0.5, float('inf')):'#08306b', 
             (0.3, 0.5):'#00441b',
             (0.1, 0.3):'#ec7014',
             (0, 0.1):'#7f0000',
             'noFrip':'#969696'}

####################################  CODE  ############################################
# Create output folders
if not os.path.exists(outplot):
    os.makedirs(outplot)
if not os.path.exists(outBigW):
    os.makedirs(outBigW)
    
# Get list of available BAM files
species = [f for f in os.listdir(inPath) if os.path.isdir(f"{inPath}/{f}")]

##################################
# Load Sample data table
##################################
df_data = pd.read_table(sample_table)
pos_control = df_data['is.control'] == "Y"
df_data["factor"] = df_data["sample_id"].str.split('_').str[1]

##################################
# Get bigwigs at given resolution
##################################
for sps_id in species:
    files = [f for f in os.listdir(f"{inPath}/{sps_id}") if f.endswith('bam')]
    for fi in files:
        if not os.path.exists(f"{outBigW}/{sps_id}__{fi[:-3]}bw"):
            print(f"\n{fi} -- ")
            cmd = f"{condaEnv}/bamCoverage --binSize {resolution} \
                --normalizeUsing CPM --exactScaling \
                -b {inPath}/{sps_id}/{fi} -of bigwig \
                -o {outBigW}/{sps_id}__{fi[:-3]}bw --numberOfProcessors {nCPU}"
            _ = runCommand(cmd)

# In case we have external controls, get the ones with full path
allControlsPath = df_data[["sps_id", "control"]].drop_duplicates()
allControlsPath = allControlsPath[allControlsPath["control"].isnull() == False]
allControlsPath = allControlsPath[allControlsPath['control'].apply(
                                        lambda x:os.path.isabs(x))]
#allControlsPath = allControlsPath.set_index('sps_id')['control'].to_dict()
for idx in allControlsPath.index:
    sps_id = allControlsPath.loc[idx, "sps_id"]
    fi = allControlsPath.loc[idx, "control"]
    basefi = os.path.basename(fi)
    if not os.path.exists(f"{outBigW}/{sps_id}__{basefi[:-3]}bw"):
        print(f"\n{basefi} -- ")
        cmd = f"{condaEnv}/bamCoverage --binSize {resolution} \
            --normalizeUsing CPM --exactScaling \
            -b {fi} -of bigwig \
            -o {outBigW}/{sps_id}__{basefi[:-3]}bw --numberOfProcessors {nCPU}"
        _ = runCommand(cmd)


##################################
# Get alignment and peak calling stats
##################################

# FRiP scores
df_stats = pd.read_table(gatheredResults)
df_stats.index = df_stats["sampleName"]

FripData = {}
for id1 in df_stats.index:
    FripData[id1] = float(df_stats.loc[id1, ["narrowFRiP", "broadFRiP"]].max())
    
# Read duplicates
dupRate = {}
for id1 in df_stats.index:
    dupRate[id1] = round(float(df_stats.loc[id1, "dupRate"].split("%")[0]), 2)
    dupRate[id1] = f"{dupRate[id1]}%"


##################################
# Get read counts 
##################################
# Very very fast
allBamReads = {}
for sps_id in species:
    files = [f for f in os.listdir(f"{inPath}/{sps_id}") if f.endswith('bam')]
    for bam in files:
        filename = f"{inPath}/{sps_id}/{bam}"
        # Count only mapped reads (3rd column)
        id1 = bam.split('.')[0]
        allBamReads[id1] = sum(int(l.split('\t')[2]) 
                                for l in pysam.idxstats(filename).rstrip('\n').split('\n'))


##################################
# Get distributions
##################################

files = [f for f in os.listdir(outBigW) if f.endswith('bw')]
frips = {}
colors = {}


for sps_id in species:
    if saveFig:
        pdf = matplotlib.backends.backend_pdf.PdfPages(f'{outplot}/{sps_id}_cpmDistrib_{resolution//1000}kb.pdf')
        
    print ('#' * 4, sps_id, '#' * 4)
    spsFiles = [f for f in files if f.startswith(f"{sps_id}__")]
    IDs = dict((s.split("__")[1].split('.')[0], s) for s in spsFiles)
    IDs_rev = dict((s, s.split("__")[1].split('.')[0]) for s in spsFiles)
    pos_chips = (df_data["sample_id"].isin(IDs.keys())) & (pos_control == False)
    chips = sorted(list(set(df_data.loc[pos_chips, "factor"])))
    
    # Get the counts for all the controls we are going to use
    allControls = list(set([c for c in df_data.loc[pos_chips, "control"] 
                        if ((c == c) & (c is not None))]))
    controlCounts = {}
    for ctrl in allControls:
        # If we're dealing with an external control
        if os.path.isabs(ctrl):
            baseCtrol = os.path.basename(ctrl)
            bw1 = pyBigWig.open(f"{outBigW}/{IDs[baseCtrol.split('.')[0]]}")
        else:
            bw1 = pyBigWig.open(f"{outBigW}/{IDs[ctrl]}")
        controlCounts[ctrl] = []
        useChroms = list(bw1.chroms().keys())
        for chrom in useChroms:
            nb = bw1.chroms()[chrom]//resolution + 1
            controlCounts[ctrl] += [bw1.intervals(chrom, i * resolution, (i+1) * resolution)[0][2]
                                for i in range(nb - 1)]        
        controlCounts[ctrl] = [c for c in controlCounts[ctrl] 
                                        if (c != 0 and c == c)]
        # Add a value of one to avoid log issues
        controlCounts[ctrl] = [c + 1 for c in controlCounts[ctrl]]
            
    for chip in chips:
        id_chip = df_data.loc[((df_data["factor"] == chip) & (pos_chips)), 
                                "sample_id"]
        bws = [IDs[c] for c in id_chip]
        counts = {}
        for bw in bws:
            # get plot color
            id1 = IDs_rev[bw]
            if id1 in FripData:
                frip = FripData[id1]
                frips[bw] = round(frip, 3)
                for fp in fripColors:
                    if fp != 'noFrip':
                        if fp[0] <= frip < fp[1]:
                            colors[bw] = fripColors[fp]
            else:
                print(f"{id1} not in FripData")
                colors[bw] = fripColors['noFrip']
                frips[bw] = 'NA'
            
            
            bw1 = pyBigWig.open(f"{outBigW}/{bw}")
            counts[bw] = []
            
            for chrom in useChroms:
                nb = bw1.chroms()[chrom]//resolution + 1
                counts[bw] += [bw1.intervals(chrom, i * resolution, (i+1) * resolution)[0][2]
                                    for i in range(nb - 1)]
            
            counts[bw] = [c for c in counts[bw] if (c != 0 and c == c)]
            counts[bw] = [c + 1 for c in counts[bw]]
        
        
        nbw = len(bws)
        # Make text figure
        areas = {}
        fig, ax = plt.subplots(1,nbw, figsize=(5*nbw , 5))
        if nbw == 1:
            id1 = IDs_rev[bws[0]]
            # Find matching control
            pos = df_data['sample_id'] == IDs_rev[bws[0]]
            controlID = df_data.loc[pos, "control"].item()
            if controlID in controlCounts:
                x0 = np.array(np.log2(controlCounts[controlID]))
            # No control specified for this file
            else:
                print(f"No control specified for {id1}")
                print(controlID)
                print()
                x0 = np.array([0, 1])
            x1 = np.array(np.log2(counts[bws[0]]))
            areas[id1] = gaussAreaRight(x0, x1, ax)
 
            title = id1
            title += f"\nmaxFRiP: {FripData[id1]}, reads: {allBamReads[id1]:,}"
            title += f"\n%Dupli:{dupRate[id1]}, AreaPass: {areas[id1]}"
            ax.set_title(title)
            # Identify labels
            ax.set_ylabel("KDE")
            ax.set_xlabel("Log2(Non zero CPM + 1)")

        else:
            for nb, bw in enumerate(bws):
                id1 = IDs_rev[bw]
                # Find matching control
                pos = df_data['sample_id'] == IDs_rev[bws[0]]
                controlID = df_data.loc[pos, "control"].item()
                if controlID in controlCounts:
                    x0 = np.array(np.log2(controlCounts[controlID]))
                # No control specified for this file
                else:
                    print(f"No control specified for {id1}")
                    print(controlID)
                    print()
                    x0 = np.array([0, 1])
                x1 = np.array(np.log2(counts[bw]))
                areas[id1] = gaussAreaRight(x0, x1, ax[nb])

                title = id1
                title += f"\nmaxFRiP: {FripData[id1]}, reads: {allBamReads[id1]:,}"
                title += f"\n%Dupli:{dupRate[id1]}, AreaPass: {areas[id1]}"
                ax[nb].set_title(title)
                # Identify labels
                if nb == 0:
                    ax[nb].set_ylabel("KDE")
                ax[nb].set_xlabel("Log2(Non zero CPM + 1)")
        fig.suptitle(f"{sps_id} -- {chip}", y = 1.05)
        plt.show()
        if saveFig:
            pdf.savefig(fig , bbox_inches='tight')
        plt.close()

    if saveFig:
        pdf.close()
