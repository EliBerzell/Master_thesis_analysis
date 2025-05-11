# Input list of pileup files and output files here!
pileuplist = ["Mutator_data/Single_nucleotide_analysis/A016-plus_CA-ATTAGCTGAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A017-WT-1-CTTGGTTGAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A018-minus_CA-1-TCGCCAGAAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A019-WT_2-GACTCGCAAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A020-plus_CA-2-GTTAGGAAAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A021-minus_CA-2-ACACAGGCAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A022-WT-3-CATGTTGAAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A023-plus_CA-3-CGACTGCTAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv",
              "Mutator_data/Single_nucleotide_analysis/A024-minus_CA-3-TTGCATGTAT-GTTCAGAGTT-Guleycan-0037-A_R1_sorted_pileup.tsv"]

namelist = ["Mutator_data/Single_nucleotide_analysis/A016-plus_CA-ATTAGCTGAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A017-WT-1-CTTGGTTGAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A018-minus_CA-1-TCGCCAGAAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A019-WT_2-GACTCGCAAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A020-plus_CA-2-GTTAGGAAAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A021-minus_CA-2-ACACAGGCAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A022-WT-3-CATGTTGAAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A023-plus_CA-3-CGACTGCTAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv",
            "Mutator_data/Single_nucleotide_analysis/A024-minus_CA-3-TTGCATGTAT-GTTCAGAGTT-Guleycan-0037-A_R1_mutations.tsv"]
            
annotation_file = 'Genomes/Mus_musculus/annotation.tsv'

# No edits after this point

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

annotation = pd.read_csv(annotation_file, sep='\t',
                         dtype = {'Gene': str, 'Start': int, 'End': int, 'Orientation' : str})

for file, name in zip(pileuplist, namelist):

    print('Creating '+name)
    
    df = pd.read_csv(file, sep="\t", names=["Ref","Pos","RefNT","Depth","Match","Qual"])

    poslist = []
    refntlist = []
    depthlist = []
    Alist = []
    Clist = []
    Glist = []
    Tlist = []
    modlist = []
    genelist = []

    for index, row in df.iterrows():
        gene = annotation[(annotation['Start']<=row['Pos']) & (annotation['End']>=row['Pos'])]['Gene']
        strand = annotation[(annotation['Start']<=row['Pos']) & (annotation['End']>=row['Pos'])]['Orientation']
        if ((row['Depth'] > 0) & (len(gene) > 0)):
            gene = gene.iloc[0]
            strand = strand.iloc[0]
            matchstring = row['Match']
            if strand == 'plus':
                stranddepth = (matchstring.count('A')) + (matchstring.count('C')) + (matchstring.count('G')) + (matchstring.count('T')) + (matchstring.count('.'))
            elif strand == 'minus':
                stranddepth = (matchstring.count('a')) + (matchstring.count('c')) + (matchstring.count('g')) + (matchstring.count('t')) + (matchstring.count(','))
            if stranddepth > 0:
                if strand == 'plus':
                    A = (matchstring.count('A'))/stranddepth
                    C = (matchstring.count('C'))/stranddepth
                    G = (matchstring.count('G'))/stranddepth
                    T = (matchstring.count('T'))/stranddepth
                elif strand == 'minus':
                    A = (matchstring.count('a'))/stranddepth
                    C = (matchstring.count('c'))/stranddepth
                    G = (matchstring.count('g'))/stranddepth
                    T = (matchstring.count('t'))/stranddepth
                modlist.append(A+C+G+T)
                if row['RefNT'] == 'A':
                    if strand == 'plus':
                        A = (matchstring.count('.'))/stranddepth
                    elif strand == 'minus':
                        A = (matchstring.count(','))/stranddepth
                if row['RefNT'] == 'C':
                    if strand == 'plus':
                        C = (matchstring.count('.'))/stranddepth
                    elif strand == 'minus':
                        C = (matchstring.count(','))/stranddepth
                if row['RefNT'] == 'G':
                    if strand == 'plus':
                        G = (matchstring.count('.'))/stranddepth
                    elif strand == 'minus':
                        G = (matchstring.count(','))/stranddepth
                if row['RefNT'] == 'T':
                    if strand == 'plus':
                        T = (matchstring.count('.'))/stranddepth
                    elif strand == 'minus':
                        T = (matchstring.count(','))/stranddepth
                poslist.append(row['Pos'])
                refntlist.append(row['RefNT'])
                depthlist.append(stranddepth)
                Alist.append(A)
                Clist.append(C)
                Glist.append(G)
                Tlist.append(T)
                genelist.append(gene)

    if len(poslist)>0:  
        mutations = pd.DataFrame(columns=['Loc', 'RefNT', 'Depth', 'A', 'C', 'G', 'T', 'ModIndex', 'Gene'])

        mutations['Loc'] = poslist
        mutations['RefNT'] = refntlist
        mutations['Depth'] = depthlist
        mutations['A'] = Alist
        mutations['C'] = Clist
        mutations['G'] = Glist
        mutations['T'] = Tlist
        mutations['ModIndex'] = modlist
        mutations['Gene'] = genelist

        mutations.to_csv(name, sep='\t')
        print(name+ ' created')
    else:
        print(file+' had no mutations')
