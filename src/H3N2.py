import os, time
from collections import defaultdict
from nextstrain.io_util import make_dir, remove_dir, tree_to_json, write_json
from nextstrain.sequences import sequence_set, num_date
from nextstrain.tree import tree
from nextstrain.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

class flu_process(object):
    """docstring for flu_process"""
    def __init__(self, fname = 'data/H3N2_gisaid_epiflu_sequence.fasta',
                 outgroup_file='source_data/H3N2_outgroup.gb', **kwargs):
        super(flu_process, self).__init__()
        self.fname = fname
        self.kwargs = kwargs
        tmp_outgroup = SeqIO.read(outgroup_file, 'genbank')
        self.outgroup = tmp_outgroup.features[0].qualifiers['strain'][0]
        genome_annotation = tmp_outgroup.features
        ref_seq = SeqIO.read(outgroup_file, 'genbank')
        self.proteins = {f.qualifiers['gene'][0]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                for f in ref_seq.features if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}

        self.time_interval = [datetime.strptime('2008-01-01', "%Y-%m-%d").date(),
                              datetime.strptime('2016-01-01', "%Y-%m-%d").date()]
        self.frequencies = defaultdict(dict)
        self.pivots = np.linspace(num_date(self.time_interval[0]),
                                  num_date(self.time_interval[1]),40)

        self.seqs = sequence_set(self.fname, reference=self.outgroup)
        self.seqs.ungap()
        self.seqs.parse({0:'strain', 1:'isolate_id', 3:'passage', 5:'date', 7:'lab', 8:"accession"}, strip='_')
        self.seqs.raw_seqs['A/Beijing/32/1992'].attributes['date']='1992-01-01'
        self.seqs.raw_seqs[self.outgroup].seq=tmp_outgroup.seq
        self.seqs.parse_date(["%Y-%m-%d"], prune=True)

    def subsample(self):
        self.seqs.raw_seqs = {k:s for k,s in self.seqs.raw_seqs.iteritems() if
                                        s.attributes['date']>=self.time_interval[0] and
                                        s.attributes['date']<self.time_interval[1]}
        self.seqs.subsample(category = lambda x:(x.attributes['date'].year,x.attributes['date'].month),
                            threshold=params.viruses_per_year)

    def align(self):
        self.seqs.align()
        self.seqs.strip_non_reference()
        self.seqs.clock_filter(n_iqd=3, plot=True, max_gaps=0.05, root_seq=self.outgroup)
        self.seqs.translate(proteins=self.proteins)

    def estimate_mutation_frequencies(self):
        time_points = [x.attributes['num_date'] for x in self.seqs.aln]
        aln_frequencies = alignment_frequencies(self.seqs.aln, time_points,
                                self.pivots, ws=len(time_points)/10, **self.kwargs)
        aln_frequencies.mutation_frequencies(min_freq=0.1)
        self.frequencies['nuc'] = aln_frequencies.frequencies
        for prot in self.seqs.translations:
            aln_frequencies = alignment_frequencies(self.seqs.translations[prot], time_points,
                                            self.pivots, ws=len(time_points)//10, **self.kwargs)
            aln_frequencies.mutation_frequencies(min_freq=0.01)
            self.frequencies[prot] = aln_frequencies.frequencies

    def estimate_tree_frequencies(self):
        tree_freqs = tree_frequencies(self.tree.tree, self.pivots,
                                      ws = self.tree.tree.count_terminals()//10, **self.kwargs)
        tree_freqs.estimate_clade_frequencies()
        self.tree_frequencies = tree_freqs.frequencies


    def build_tree(self):
        self.tree = tree(aln=self.seqs.aln) #, proteins = self.seqs.proteins)

        self.tree.build(root='oldest')
        self.tree.ancestral()
        self.tree.timetree(Tc=0.1)
        self.tree.add_translations()
        self.tree.refine()
        self.tree.layout()

    def export(self, prefix='web/data/'):
        def process_freqs(freq):
            return [round(x,4) for x in freq]
        self.seqs.export_diversity(prefix+'entropy.json')
        self.tree.export(path=prefix, extra_attr = ["subtype", "country", "region", "nuc_muts", "aa_muts"])

        freq_json = {}
        for gene, tmp_freqs in self.frequencies.iteritems():
            for mut, freq in tmp_freqs.iteritems():
                freq_json['_'.join([gene, str(mut)[0], mut[1]])] = process_freqs(freq)
        for clade, freq in self.tree_frequencies.iteritems():
            freq_json['clade_'+str(clade)] = process_freqs(freq)
        write_json(freq_json, prefix+'frequencies.json', indent=None)



if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-g', '--gene', type = str, default='pol', help='The HIV gene to process')
    parser.add_argument('-s', '--subtype', type = str, default='B', help='The HIV subtype to filter')
    parser.add_argument('-v', '--viruses_per_year', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')

    params = parser.parse_args()

    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9)
    flu.subsample()
    flu.align()
    flu.estimate_mutation_frequencies()
    flu.build_tree()
    flu.estimate_tree_frequencies()
    flu.export()
