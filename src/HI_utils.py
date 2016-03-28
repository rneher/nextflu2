import numpy as np
import time, os, gzip
from collections import defaultdict
from nextstrain.io_util import myopen
from matplotlib import pyplot as plt
from itertools import izip
from H3N2 import fix_name
import pandas as pd


######################################################################################
####  utility functions for reading and writing HI Tables, plotting trees, etc
######################################################################################


def HI_fix_name(name):
    if name.split() == ["NIB-85", "(A/Almaty/2958/2013)"]:
        tmp_name = fix_name("A/Almaty/2958/2013")
    elif name.split() == ["A/Texas/50/2012","(6&7)"]:
        tmp_name = fix_name("A/Texas/50/2012")
    else:
        tmp_name = fix_name(name)
    return tmp_name.upper().lstrip('*')


def plot_tree(tree):
    from Bio import Phylo
    from tree_util import to_Biopython, color_BioTree_by_attribute
    btree = to_Biopython(tree)
    color_BioTree_by_attribute(btree,"cHI", transform = lambda x:x)
    Phylo.draw(btree, label_func = lambda  x: 'X' if x.serum else '' if x.HI_info else '',
        show_confidence= False) #, branch_labels = lambda x:x.mutations)

def plot_dHI_distribution(tree):
    plt.figure()
    mut_dHI = sorted([(c.dHI, c.n_aa_muts) for c in tree.find_clades(order='postorder') if c.HI_info],
                    reverse=True, key=lambda x:x[0])
    thres = [0,1,2,3,100]
    for lower, upper in zip(thres[:-1], thres[1:]):
        tmp = [x[0] for x in mut_dHI if len(x[1])>=lower and len(x[1])<upper]
        plt.plot(sorted(tmp), np.linspace(1,0, len(tmp)), label='#aa='+str(lower)+', total '+str(len(tmp)))
    plt.legend(loc=1)
    plt.show()

min_titer = 10.0

def titer_to_number(val):
    try:
        if '<' in val:
            return np.nan
        if len(val.split())>1:
            return float(val.split()[0])
        else:
            return float(val)
    except:
        #print "Bad HI measurement:", val
        return np.nan

def parse_HI_matrix(fname):
    from string import strip
    import csv
    name_abbrev = {'HK':"HONGKONG", 'SWITZ':"SWITZERLAND", 'VIC':"VICTORIA", 'STOCK':"STOCKHOLM",
                    'STHAFR':"SOUTHAFRICA", 'SAFRICA':"SOUTHAFRICA", "ENG":"ENGLAND", "NIB-85":"A/ALMATY/2958/2013", 'NOR':'NORWAY',
                    'NTHCAROL':"NORTHCAROLINA",'ALA':"ALABAMA", 'NY':"NEWYORK", "GLAS":"GLASGOW", "AL":"ALABAMA",
                    "NETH":"NETHERLANDS", "FIN":"FINLAND", "BRIS":"BRISBANE", "MARY":"MARYLAND",
                    "ST.P'BURG":"ST.PETERSBURG", 'CAL':'CALIFORNIA', 'AUCK':'AUCKLAND', "C'CHURCH":'CHRISTCHURCH',
                    'CHCH':'CHRISTCHURCH', 'ASTR':'ASTRAKHAN', 'ASTRAK':'ASTRAKHAN', 'ST.P':"ST.PETERSBURG",'ST P':"ST.PETERSBURG",'STP':"ST.PETERSBURG",
                    'JHB':'JOHANNESBURG', 'FOR':'FORMOSA','MAL':'MALAYSIA', 'STHAUS':'SOUTHAUSTRALIA',
                    'FL':'FLORIDA', 'MASS':'MASSACHUSETTS','NOVO':'NOVOSIBIRSK','WIS':'WISCONSIN','BANG':'BANGLADESH','EG':'EGYPT'  }
    src_id = fname.split('/')[-1]
    print fname
    with myopen(fname) as infile:
        csv_reader = csv.reader(infile)

        # parse sera
        row1 = csv_reader.next()
        row2 = csv_reader.next()
        row3 = csv_reader.next()
        ref_sera = [[HI_fix_name(e1+'/'+e2), e3.replace(' ','')] for e1,e2,e3 in zip(row1, row2, row3)[4:]]
        for ri in xrange(len(ref_sera)):
            abbr = ref_sera[ri][0].split('/')[1].rstrip('01234566789')
            if abbr in name_abbrev:
                ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0].replace(abbr, name_abbrev[abbr]))
            else:
                ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0])
            # strip numbers
            tmp = ref_sera[ri][0].split('/')
            ref_sera[ri][0] = '/'.join([tmp[0], tmp[1].rstrip('0123456789')]+tmp[2:])
            try:
                y = int(ref_sera[ri][0].split('/')[-1])
                if y<100:
                    if y<20:
                        ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(2000+y)
                    else:
                        ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(1900+y)
            except:
                print ref_sera[ri]

        fields = ['source','ref/test', 'genetic group', 'collection date', 'passage history']+map(tuple, ref_sera)
        #print fields
        for row in csv_reader: # advance until the reference virus
            if row[0].startswith('REFERENCE'):
                break

        ref_strains = []
        ref_matrix = []
        for row in csv_reader:
            if row[0].startswith('TEST'):
                break
            else: # load matrices until the test virus section starts
                ref_strains.append(HI_fix_name(row[0].strip()))
                ref_matrix.append([src_id,'ref']+map(strip, row[1:4])+map(titer_to_number, row[4:]))

        test_strains = []
        test_matrix = []
        for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
            if not (row[0].startswith('A/') or row[0].startswith('B/')):
                break
            else:
                test_strains.append(HI_fix_name(row[0].strip()))
                test_matrix.append([src_id,'test']+map(strip,row[1:4])+map(titer_to_number, row[4:]))

        print len(ref_sera), ref_sera
        print len(ref_strains), len(test_strains)
        HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)

        return HI_table


def read_tables(flutype = 'H3N2'):
    import glob
    from itertools import product
    flist = glob.glob('../../HI_titers/'+flutype+'_tables/NIMR*csv')
    all_names = set()
    all_measurements = defaultdict(list)
    HI_matrices = pd.DataFrame()
    for fname in flist:
        tmp = parse_HI_matrix(fname)
        HI_matrices = HI_matrices.append(tmp)
    return HI_matrices

def read_trevor_table(flutype):
    trevor_table = 'data/'+flutype+'_HI.tsv'
    import csv
    measurements = []
    sera = set()
    strains = set()
    if os.path.isfile(trevor_table):
        with myopen(trevor_table) as infile:
            table_reader = csv.reader(infile, delimiter="\t")
            header = table_reader.next()
            for row in table_reader:
                val = titer_to_number(row[6])
                if not np.isnan(val):
                    strains.add(HI_fix_name(row[1]))
                    serum = (HI_fix_name(row[4]), row[3])
                    src_id = row[-1]
                    sera.add(serum)
                    measurements.append([HI_fix_name(row[1]), serum, src_id, val])
    else:
        print trevor_table, "not found"
    print "trevor total:", len(measurements), "measurements"
    return strains, sera, pd.DataFrame(measurements)


def table_to_flat(HI_table):
    flat_measurements = list()
    for ref_serum in HI_table.columns[5:]:
        sub_set_vals = HI_table[ref_serum][~np.isnan(HI_table[ref_serum])]
        sub_set_source = HI_table['source'][~np.isnan(HI_table[ref_serum])]
        for virus, val, src_id in izip(sub_set_vals.index, sub_set_vals, sub_set_source):
            flat_measurements.append([virus, ref_serum, src_id, val])
    print "NIMR total:", len(flat_measurements), "measurements"
    return pd.DataFrame(flat_measurements)

def get_all_titers_flat(flutype='H3N2'):
    HI_titers = read_tables(flutype)
    HI_titers_flat = table_to_flat(HI_titers)
    HI_trevor = read_trevor_table(flutype)[2]
    HI_titers_flat = pd.concat((HI_titers_flat, HI_trevor))
    return HI_titers_flat


def write_strains_with_HI_and_sequence(flutype='H3N2'):
    HI_titers = read_tables(flutype)
    HI_trevor = read_trevor_table(flutype)
    HI_strains = set(HI_titers.index)
    HI_strains.update(HI_trevor[0])
    from Bio import SeqIO
    good_strains = set()
    with myopen("data/"+flutype+"_strains_with_HI.fasta", 'w') as outfile, \
         myopen("data/"+flutype+"_gisaid_epiflu_sequence.fasta", 'r') as infile:
        for seq_rec in SeqIO.parse(infile, 'fasta'):
            tmp_name = seq_rec.description.split('|')[0].strip()
            reduced_name = HI_fix_name(tmp_name)
            if reduced_name in HI_strains and (reduced_name not in good_strains):
                SeqIO.write(seq_rec, outfile,'fasta')
                good_strains.add(reduced_name)

    titer_count = defaultdict(int)
    measurements = get_all_titers_flat(flutype)
    for ii, rec in measurements.iterrows():
        test, ref, src_id, val = rec
        titer_count[test]+=1

    with myopen("data/"+flutype+"_HI_strains.txt", 'w') as HI_strain_outfile:
        for strain, count in sorted(titer_count.items(), key=lambda x:x[1], reverse=True):
            HI_strain_outfile.write(strain + '\t'+str(count)+'\n')
            if fix_name(strain)!=strain:
                HI_strain_outfile.write(fix_name(strain) + '\t'+str(count)+'\n')


def write_flat_HI_titers(flutype = 'H3N2', fname = None):
    measurements = get_all_titers_flat(flutype)
    with myopen('data/'+flutype+'_HI_strains.txt') as infile:
        strains = [HI_fix_name(line.strip().split('\t')[0]).upper() for line in infile]
    if fname is None:
        fname = 'data/'+flutype+'_HI_titers.txt'
    written = 0
    skipped = 0
    with myopen(fname, 'w') as outfile:
        for ii, rec in measurements.iterrows():
            test, ref, src_id, val = rec
            if HI_fix_name(test).upper() in strains and HI_fix_name(rec[1][0]).upper() in strains:
                outfile.write('\t'.join(map(str, [test, ref[0], ref[1], src_id, val]))+'\n')
                written+=1
            else:
                skipped+=1
    print "written",written,"records"
    print "skipped",skipped,"records"

def export_ref_strains(myflu):
    strains = []
    for r in myflu.ref_strains:
        tmp = myflu.sequence_lookup[myflu.node_lookup[r].strain]
        strains.append({'seq': str(tmp.seq), 'date': tmp.date, 'strain': tmp.strain, 'region': tmp.region, 'country':tmp.country})
    from json import dump as jdump
    with open('source-data/'+ myflu.virus_type+'_ref_strains.json', 'w') as ofile:
        jdump(strains, ofile, indent=2)
