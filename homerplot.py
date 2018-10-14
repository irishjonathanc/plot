#!/usr/bin/env python

# Plotting module for homertools ChIPseq quality control output
# Jonathan Irish, June 2014

import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

def plot_NucFreq(filename, title):
    """Plot mean nucleotide frequency at each base position"""

    input_freq = pd.read_csv("{0}".format(filename), sep="\t")
    a = input_freq["A"]
    c = input_freq["C"]
    g = input_freq["G"]
    t = input_freq["T"]
    offset = input_freq["Offset"]
    matplotlib.rcParams.update({'font.size' : 14})
    freq_fig, axes = plt.subplots(figsize=(8,8))
    axes.plot(offset, a, label="A", lw=2)
    axes.plot(offset, c, label="C", lw=2)
    axes.plot(offset, g, label="G", lw=2)
    axes.plot(offset, t, label="T", lw=2)
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')
    axes.yaxis.tick_left()
    axes.xaxis.tick_bottom()
    axes.set_xlabel("distance from 5' end of reads")
    axes.set_ylabel('Nucleotide frequency')
    axes.set_title('Nucleotide frequency for {0} Sample'.format(title), fontsize=16)
    axes.legend(loc=1, fontsize=12, frameon=False)
    freq_fig.savefig("{0}_NucFreqPlot.pdf".format(title))


def plot_TagCorr(filename, title):
    """Plot tag correlation"""

    df=pd.read_csv(filename, sep="\t", skiprows=1, header=None, names=['distance', 'plus_strand', 'minus_strand'])
    x1 = df["distance"]
    y1 = df["plus_strand"]
    z1 = df["minus_strand"]
    matplotlib.rcParams.update({'font.size' : 14})
    tag_fig, axes = plt.subplots(figsize=(8,5))
    axes.plot(x1, y1, label="plus_strand")
    axes.plot(x1, z1, label="minus_strand")
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')
    axes.yaxis.tick_left()
    axes.xaxis.tick_bottom()
    axes.set_xlabel('distance')
    axes.set_ylabel('reads')
    axes.set_title('Tag Autocorrelation for {}'.format(title), fontsize=14)
    axes.legend(loc=1, fontsize=12, frameon=False);
    tag_fig.savefig("{0}_tagCorrPlot.pdf".format(title))

def plot_TagDup(filename, title):
    """Plot tag duplication levels"""

    tagDist = pd.read_csv(filename, sep = '\t') 
    i = tagDist["clones"]
    j = tagDist["frequency"]
    dup_fig, ax = plt.subplots(figsize=(8,5))
    ax.bar(i, j, align = "center", width = 0.5, alpha = 0.8)
    ax.set_xlim([0,10])
    ax.set_xticks(range(11))
    ax.set_xlabel('Reads per position')
    ax.set_ylabel('Fraction of total reads')
    ax.set_title('Tag Duplication Frequency for {}'.format(title), fontsize = 18)
    dup_fig.savefig("{0}_DupPlot.pdf".format(title))

def plot_genomeGC_dist(genome_file, sample_file, title):

    genome = pd.read_csv(genome_file, sep = '\t')
    sample = pd.read_csv(sample_file, sep = '\t')
    genome_gc = genome["GC%"]
    genome_fraction = genome["Fraction"]
    sample_gc = sample["GC%"]
    sample_fraction = sample["Fraction"]
    matplotlib.rcParams.update({'font.size' : 14})
    gc_fig, axes = plt.subplots(figsize=(10,8))
    axes.plot(genome_gc, genome_fraction, label = "genome", lw=4)
    axes.plot(sample_gc, sample_fraction, label = "{}".format(title), lw=4)
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')
    axes.yaxis.tick_left()
    axes.set_xlabel("Genome %GC content")
    axes.set_ylabel('Fraction of Genome')
    axes.set_title('Genome %GC content for {0}'.format(title), fontsize=16)
    axes.legend(loc=1, fontsize=12, frameon=False)
    gc_fig.savefig("{0}_Genome_%GC_Plot.pdf".format(title))

if __name__ == "__main__":

    # parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("tagfile_in", help = "Enter a tag correlation filename")
    parser.add_argument("dupfile_in", help = "Enter a count distribution filename")
    parser.add_argument("gcfreq_in", help = "Enter a tagFreq filename")
    parser.add_argument("genomeGC_in", help = "Enter a genomeGC filename")
    parser.add_argument("sample_name", help = "Enter a sample name for output files and\
        plot titles")
    args = parser.parse_args()
    print "Input files: {} {} {}".format(args.tagfile_in, args.dupfile_in, \
      args.gcfreq_in, args.genomeGC_in) 
    print "Sample name is {}".format(args.sample_name)

    # generate plots and save to .pdf
    plot_NucFreq(args.gcfreq_in, args.sample_name)
    plot_TagCorr(args.tagfile_in, args.sample_name)
    plot_TagDup(args.dupfile_in, args.sample_name)
    plot_genomeGC_dist(args.genomeGC_in, args.sample_name)
