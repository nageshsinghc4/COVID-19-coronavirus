#task to perform:
# 1. Protein Analysis
# 2. Comparing Human Coronavirus RNA with MERS and SARS


#pip install seaborn
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
pd.plotting.register_matplotlib_converters()
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns

import os
"""
Some facts regarding the virus itself: COVID-19

It is of the (+)ssRNA classification of viruses, which means it is a single stranded virus that can be directly translated into protein.
The actual virus is called SARS-CoV-2, Covid-19 is the name for the respiratory disease it causes (I found this interesting)
"""

#pip install biopython
from Bio import SeqIO
for sequence in SeqIO.parse('/Users/nageshsinghchauhan/Downloads/coronavirus-genome-sequence/MN908947.fna', "fasta"):
    print(sequence.id)
    print(sequence.seq)
    print(len(sequence),'nucliotides')

# Loading Complementary DNA Sequence into an alignable file
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
DNAsequence = SeqIO.read('/Users/nageshsinghchauhan/Downloads/coronavirus-genome-sequence/MN908947.fna', "fasta")

DNA = DNAsequence.seq


#Convert DNA into mRNA Sequence 
mRNA = DNA.transcribe() #Transcribe a DNA sequence into RNA.
print(mRNA)
print(len(mRNA))

"""
The difference between the complementary DNA and the mRNA is just that the bases T (for Thymine) is replaced with U (for Uracil).
"""
# Obtain Amino Acid Sequence from mRNA
Amino_Acid = mRNA.translate(table=1, cds=False)
print('Amino Acid : ', Amino_Acid)
print("Length of Protein : ",len(Amino_Acid))
print("Length of Original mRNA : ",len(mRNA))

"""
Codons
Cells decode mRNAs by reading their nucleotides in groups of three, called codons. Here are some features of codons:
Most codons specify an amino acid
Three "stop" codons mark the end of a protein
One "start" codon, AUG, marks the beginning of a protein and also encodes the amino acid methionine
"""
from Bio.Data import CodonTable
print(CodonTable.unambiguous_rna_by_name['Standard'])

"""
Let's now identify all the polypeptides so basically separating at the stop codon, marked by * . Then let's remove any sequence less than 20 amino acids long, as this is the smallest known functional protein (if curious). Note: In humans the smallest known functional protien is 44 amino acids long.
"""
#Identify all the Proteins (chains of amino acids)
Proteins = Amino_Acid.split('*') # * is translated stop codon
for i in Proteins[:]:
    if len(i) < 20:
        Proteins.remove(i)

# 1. Protein Analysis With The Protparam Module In Biopython
from __future__ import division
poi_list = []
MW_list = []
from Bio.SeqUtils import ProtParam
for record in Proteins[:]: 
    print("\n")
    X = ProtParam.ProteinAnalysis(str(record))
    POI = X.count_amino_acids()
    poi_list.append(POI)
    MW = X.molecular_weight()
    MW_list.append(MW)
    print("Protein of Interest = ", POI) 
    print("Amino acids percent = ", str(X.get_amino_acids_percent())) 
    print("Molecular weight = ", MW)
    print("Aromaticity = ", X.aromaticity()) 
    print("Flexibility = ", X.flexibility()) 
    print("Isoelectric point = ", X.isoelectric_point()) 
    print("Secondary structure fraction = ", X.secondary_structure_fraction())


MoW = pd.DataFrame(data = MW_list,columns = ["Molecular Weights"] )
MoW.head()

#plot POI
poi_list = poi_list[48]
plt.figure(figsize=(10,6));
plt.bar(poi_list.keys(), list(poi_list.values()), align='center')

# Plot lengths
plt.figure(figsize=(20,5))
plt.subplot(111)
plt.hist(functional_proteins['length'])
plt.title('Length of proteins -- histogram')
# Remove the extremes
plt.figure(figsize=(20,5))
wo = functional_proteins.loc[functional_proteins['length'] < 60]
plt.subplot(121)
plt.hist(wo['length'])
plt.title('Lenght of proteins (where < 60)')

wo = functional_proteins.loc[functional_proteins['length'] > 1000]
plt.subplot(122)
plt.hist(wo['length'])
plt.title('Length of proteins (where > 1000)')

# See what's about that huge protein
large_prot = functional_proteins.loc[functional_proteins['length'] > 2700]
l = large_prot['sequence'].tolist()[0]
print('Sequence sample:', '...',l[1000:1150],'...')

#2. Comparing Human Coronavirus RNA¶

from Bio import pairwise2
# Define sequences to be aligned
SARS = SeqIO.read("/Users/nageshsinghchauhan/Downloads/coronavirus-genome-sequence/sars.fasta", "fasta")
MERS = SeqIO.read("/Users/nageshsinghchauhan/Downloads/coronavirus-genome-sequence/mers.fasta", "fasta")
COV2 = SeqIO.read("/Users/nageshsinghchauhan/Downloads/coronavirus-genome-sequence/cov2.fasta", "fasta")

print('Sequence Lengths:')
print('SARS:', len(SARS.seq))
print('COV2:', len(COV2.seq))
print('MERS:', len(MERS.seq))

#Visualize DNA sequence using Squiggle #run in terminal
Squiggle cov2.fasta sars.fasta mers.fasta --method=gates --separate

# Alignments using pairwise2 alghoritm
SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)
MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)
MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)

# Plot the data
X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
plt.title('Sequence identity (%)')
plt.bar(X,Y)
