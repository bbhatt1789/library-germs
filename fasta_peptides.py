from Bio import SeqIO, Entrez
import pandas as pd
import numpy as np
from pyteomics import parser
from pyteomics.mass import mass
from collections import defaultdict

date_ = time.strftime('%m/%d/%Y %H:%M:%S %p')
#subprocess.call('RefProtDB download 9606 10090')

mass_dict = dict(mass.std_aa_mass)
mass_dict['X'] = np.mean([x for x in mass_dict.values()])
mass_dict['B'] = np.mean([mass_dict['D'], mass_dict['N']])
mass_dict['J'] = np.mean([mass_dict['I'], mass_dict['L']])
mass_dict['Z'] = np.mean([mass_dict['E'], mass_dict['Q']])
mass_dict['U'] = 167.057 # SELENOCYSTEINE
mass_dict['O'] = 255.318 # pyrrolysine

def input_fasta(fas = None): #read fasta
    seq = SeqIO.parse(fas, 'fasta')
    headers = []
    sequence = []
    refac = []
    for ix in seq:
        if   ix.description.count('gi') == 1: #select unique protein
            
            headers.append(ix.description.split('|'))
            sequence.append(str(ix.seq))

    for x in headers:
        refac.append(str(x[3]).split('_'))
    ref = pd.DataFrame(refac)   
    fas = pd.DataFrame(headers, sequence).reset_index()
    fasta = pd.concat([ref,fas], axis = 1) #merge on 
    fasta.columns = ['np/xp', 'np/xp_no','sequence','gi','protein_gi','ref','ref_acc','description']
    return fasta

def digest_protein(protein = None):
    digest = parser.cleave(protein, '[KR]', missed_cleavages = 3, min_length = 6)              
    return digest
    

    
def digest(ref = None):
    pep_ = pd.DataFrame()
    for ix in splitDataFrameIntoSmaller(ref):
        ix['peptide'] = ix.apply(lambda x : ','.join(digest_protein(x['sequence'])), axis =1)
        pep_ = pep_.append(ix.peptide.str.split(",").apply(pd.Series, 1).stack().reset_index())  
    return pep_
    
    
def merge_pep(pep_ = None, ref = None):    
    ref2pep = pd.merge(pep_, ref, left_on='level_0', right_index=True)
    ref2pep.rename(columns= {0 :'pep_Sequence'}, inplace = True)
    return ref2pep
    
def proteinlist(ref2pep_ = None):
    plist_tgs = ref2pep_.groupby(['pep_Sequence'])['protein_gi'].agg(lambda x: ','.join(set(x))).reset_index()
    
    plist_tgs['pep_ProteinCount'] = plist_tgs['protein_gi'].apply(lambda x: len(x.split(',')))
    plist_tgs.rename(columns = {'protein_gi': 'pep_ProteinList'}, inplace = True)
    return plist_tgs.drop_duplicates()
    
def main():
    refseq = input_fasta(fas = 'ribosomal_Bacterial_2016_ribosome.fasta')   
    pep = digest(refseq)
    ref2pep = merge_pep(pep_ = pep, ref = refseq )
    peptide = proteinlist(ref2pep_ = ref2pep)
    all_pep = pd.merge(peptide, ref2pep, on = 'pep_Sequence', how = 'left')
    all_pep.drop(columns = ['level_0','level_1'], inplace = True)
    return all_pep.drop_duplicates()

peptides = main()    
peptides.loc[peptides['pep_ProteinCount'] == 1].to_csv('unique_ribosomal_bacterial_proteins.tab', sep = '\t', index = False)