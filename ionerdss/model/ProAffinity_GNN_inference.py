# Requirements: NOTE: Change the torch version according to your CUDA version
# python >=3.8
# Step 1: Use pip to install torch==2.2.2 torchvision==0.17.2 torchaudio==2.2.2:
# go to https://docs.pytorch.org/get-started/previous-versions/, search for version 2.2.2 and find your system and GPU version
# For example, if system is Linux and GPU version is CUDA 12.1, run this:
# pip install torch==2.2.2 torchvision==0.17.2 torchaudio==2.2.2 --index-url https://download.pytorch.org/whl/cu121
# Step 2:
# pip install torch_geometric==2.3.0 transformers==4.38
#!/usr/bin/env python

# Suppress known benign warnings from dependencies
import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='huggingface_hub.file_download')
warnings.filterwarnings('ignore', message='.*Some weights of EsmModel were not initialized.*')
warnings.filterwarnings('ignore', message='.*You should probably TRAIN this model.*')

# Also suppress transformers warnings at import time
import os
os.environ['TRANSFORMERS_NO_ADVISORY_WARNINGS'] = '1'

# Suppress all transformers warnings
import logging
logging.getLogger("transformers").setLevel(logging.ERROR)

import copy
import math
import torch
import pickle
import os
import re
from itertools import chain

from torch_geometric.nn.models import AttentiveFP
from torch.optim import Adam
from torch.nn import MSELoss
from torch.nn import L1Loss
import torch.nn.functional as F

from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
import numpy as np

np.random.seed(0)
torch.manual_seed(0)


# For parsing PDB files
aminoacid_abbr = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'PHE': 'F', 'TRP': 'W', 'TYR': 'Y', 'ASP': 'D', 'ASN': 'N', 'GLU': 'E', 'LYS': 'K', 'GLN': 'Q', 'MET': 'M', 'SER': 'S', 'THR': 'T', 'CYS': 'C', 'PRO': 'P', 'HIS': 'H', 'ARG': 'R', 'UNK': 'X'}

def infer_atom_type(atom_name):
    """Infer atom type based on atom name for PDB files."""
    atom_name = atom_name.strip().upper()
    if atom_name.startswith('C'):
        return 'C'
    elif atom_name.startswith('O'):
        return 'OA'
    elif atom_name.startswith('N'):
        return 'N'
    elif atom_name.startswith('S'):
        return 'SA'
    elif atom_name.startswith('H'):
        return 'HD'
    else:
        return 'A'  # fallback or generic type

def only_letters(s):
    return re.sub('[^a-zA-Z]', '', s)

def get_chainlist_from_indexfile(chainindex):
    pdb_list = []
    chain_list = []
    with open(chainindex, 'r') as f:
        lines = f.readlines()  
        for line in lines:
            pdb = line.split('\t')[0].strip()
            pdb_list.append(pdb)
            chains = line.split('\t')[1].strip()
            chains = chains.split(';')
            while '' in chains:
                chains.remove('')

            for i, c in enumerate(chains):
                chains[i] = only_letters(c).strip()
            chain_list.append(chains)
            
    return (pdb_list, chain_list)

def find_first_numeric_part(s):
    match = re.search(r'\d+', s)
    return match.group(0) if match else None

def check_res_number(residuelist):
    char_index = []
    for i, res in enumerate(residuelist):
        if res['number'].isnumeric() == False:
            char_index.append(i)
    if len(char_index) == 0:
        return
        
    list_len = len(residuelist) 
    if(char_index[-1] == list_len - 1):
        last_num = find_first_numeric_part(residuelist[-1]['number'])
        residuelist[-1]['number'] = last_num

    for i in reversed(char_index):
        if i == list_len - 1:
            continue

        residuelist[i]['number'] = residuelist[i+1]['number']
        for j in range(i+1, len(residuelist)):
            residuelist[j]['number'] = int(residuelist[j]['number']) + 1

def get_residue_list_from_file(filename, chain_list):
    residue_chain_list = []
    try: 
        with open(filename, 'r') as f:
            lines = f.readlines()

            for chain in chain_list:
                previous_res = -1
                previous_res_type = ''
                atomlist = []
                atom = {}       
                residue = {}
                current_res_chain = ''
                residuelist = []

                for line in lines:
                    atomline = line.split()
                    if atomline[0].strip() != 'ATOM' and atomline[0].strip() != 'TER':
                        continue

                    if line[21] != chain:
                        continue

                    if atomline[0] == 'TER':
                        residue['type'] = aminoacid_abbr.get(previous_res_type, 'X')
                        residue['number'] = previous_res
                        residue['atoms'] = copy.deepcopy(atomlist)
                        residue['chain'] = current_res_chain

                        has_CA = any(atomCA['type'] == 'CA' for atomCA in residue['atoms'])
                        if has_CA:
                            residuelist.append(copy.deepcopy(residue))
                        else:
                            residue.clear()

                        previous_res = -1
                        continue

                    atom['pdbqt_type'] = infer_atom_type(line[12:16].strip())
                    atom['type'] = line[12:16].strip()
                    atom['x'] = line[30:38].strip()
                    atom['y'] = line[38:46].strip()
                    atom['z'] = line[46:54].strip()

                    current_res = line[22:27].strip()
                    current_res_type = line[17:21].strip()
                    current_res_chain = line[21]

                    if current_res != previous_res and previous_res != -1:
                        residue['type'] = aminoacid_abbr.get(previous_res_type, 'X')
                        residue['number'] = previous_res
                        residue['atoms'] = copy.deepcopy(atomlist)
                        residue['chain'] = current_res_chain

                        has_CA = any(atomCA['type'] == 'CA' for atomCA in residue['atoms'])
                        if has_CA:
                            residuelist.append(copy.deepcopy(residue))
                        else:
                            residue.clear()

                        atomlist.clear()
                        atomlist.append(copy.deepcopy(atom))
                    else:
                        atomlist.append(copy.deepcopy(atom))

                    previous_res = current_res
                    previous_res_type = current_res_type

                check_res_number(residuelist)
                residue_chain_list.append(copy.deepcopy(residuelist))

        return residue_chain_list
    except Exception as e:
        print(e)
        return residue_chain_list

# Given 2 parts of protein-protein complex, return the residue list of the complex
# Should be completed 2 residue lists

def get_interaction_residue_pair_new(dis_thred, reslistA, reslistB):
    res_pair = []
    proteinA = copy.deepcopy(reslistA)
    proteinB = copy.deepcopy(reslistB)
    
    for i, resA in enumerate(proteinA):     
        find_match_res_pair = 0
        ca1 = list(filter(lambda x: x['type'] == 'CA', resA['atoms']))[0]
        ca1_x = float(ca1['x'])
        ca1_y = float(ca1['y'])
        ca1_z = float(ca1['z'])
        
        for j, resB in enumerate(proteinB):
            find_match_res_pair = 0
            ca2 = list(filter(lambda x: x['type'] == 'CA', resB['atoms']))[0]
            ca2_x = float(ca2['x'])
            ca2_y = float(ca2['y'])
            ca2_z = float(ca2['z'])  
        
            for atomA in resA['atoms']:
                Ax = float(atomA['x'])
                Ay = float(atomA['y'])
                Az = float(atomA['z'])        

                for atomB in resB['atoms']:
                    Bx = float(atomB['x'])
                    By = float(atomB['y'])
                    Bz = float(atomB['z'])

                    distance = math.sqrt((Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2)
          
                    if distance <= dis_thred:
                        c_alpha_dist = math.sqrt((ca1_x-ca2_x)**2 + (ca1_y-ca2_y)**2 + (ca1_z-ca2_z)**2)
                        res_pair.append((resA, resB, c_alpha_dist))
                        find_match_res_pair = 1
                        break

                if find_match_res_pair == 1:
                    break

    return res_pair

def get_interaction_residue_pair_new_indi(dis_thred, reslistA, reslistB):
    res_pair = []
    proteinA = copy.deepcopy(reslistA)
    proteinB = copy.deepcopy(reslistB)
    
    for i, resA in enumerate(proteinA):     
        find_match_res_pair = 0
        ca1 = list(filter(lambda x: x['type'] == 'CA', resA['atoms']))[0]
        ca1_x = float(ca1['x'])
        ca1_y = float(ca1['y'])
        ca1_z = float(ca1['z'])
        res1_num = resA['number']
        
        for j, resB in enumerate(proteinB):
            find_match_res_pair = 0
            ca2 = list(filter(lambda x: x['type'] == 'CA', resB['atoms']))[0]
            ca2_x = float(ca2['x'])
            ca2_y = float(ca2['y'])
            ca2_z = float(ca2['z'])  
            res2_num = resB['number']
            if res1_num == res2_num:
                continue
        
            for atomA in resA['atoms']:
                Ax = float(atomA['x'])
                Ay = float(atomA['y'])
                Az = float(atomA['z'])        

                for atomB in resB['atoms']:
                    Bx = float(atomB['x'])
                    By = float(atomB['y'])
                    Bz = float(atomB['z'])

                    distance = math.sqrt((Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2)
          
                    if distance <= dis_thred:
     
                        c_alpha_dist = math.sqrt((ca1_x-ca2_x)**2 + (ca1_y-ca2_y)**2 + (ca1_z-ca2_z)**2)
                        res_pair.append((resA, resB, c_alpha_dist))
                        find_match_res_pair = 1
                        break

                if find_match_res_pair == 1:
                    break

    return res_pair

# get the edge index of the graph from the interaction pairs
# only 2 parts in the inter_pairs
# len_pro1 is the length of the first part of the protein-protein complex

def get_edge_index(inter_pairs, len_protein1_fasta):
    source_list = []
    des_list = []

    for pair in inter_pairs:
        source_res = int(pair[0]['number'])
        des_res = int(pair[1]['number']) + len_protein1_fasta
        source_list.append(source_res)
        des_list.append(des_res)

    source = torch.unsqueeze(torch.tensor(source_list), 0)
    des = torch.unsqueeze(torch.tensor(des_list), 0)

    bi_direct_source = torch.cat((source, des), 1)
    bi_direct_edge = torch.cat((des, source), 1)
    edge_index = torch.squeeze(torch.stack((bi_direct_source, bi_direct_edge)))
    
    return edge_index

def get_edge_index_indi(inter_pairs):
    source_list = []
    des_list = []

    for pair in inter_pairs:
        source_res = int(pair[0]['number'])
        des_res = int(pair[1]['number'])
        if source_res != des_res:
            source_list.append(source_res)
            des_list.append(des_res)

    source = torch.unsqueeze(torch.tensor(source_list), 0)
    des = torch.unsqueeze(torch.tensor(des_list), 0)

    edge_index = torch.squeeze(torch.stack((source, des)))
    return edge_index

def adjust_residuelist_num(residuelist):
    index = 0
    for reslist in residuelist:
        for res in reslist:
            res['number'] = index
            index += 1
            
# Helpers
def get_fasta_seq(pdb_tuple):
    fastaA = pdb_tuple[0]
    fastaB = pdb_tuple[1]
    return fastaA, fastaB

def get_distance(x1, y1, z1, x2, y2, z2):
    distance = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return distance

# Models
class AttentiveFPModel(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout):
        super(AttentiveFPModel, self).__init__()
        self.model = AttentiveFP(in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout)

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        return self.model(x, edge_index, edge_attr, batch)

class GraphNetwork(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout, linear_out1, linear_out2):
        super(GraphNetwork, self).__init__()
        self.graph1 = AttentiveFPModel(in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout)
        self.graph2 = AttentiveFPModel(in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout)
        self.graph3 = AttentiveFPModel(in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout)
        
        self.fc1 = torch.nn.Linear(out_channels * 3, linear_out1)
        self.fc2 = torch.nn.Linear(linear_out1, linear_out2)

    def forward(self, inter_data, intra_data1, intra_data2):

        inter_graph = self.graph1(inter_data)
        intra_graph1 = self.graph2(intra_data1)
        intra_graph2 = self.graph3(intra_data2)

        x = torch.cat([inter_graph, intra_graph1, intra_graph2], dim=1)
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x
        
def run_proaffinity_inference(
        pdbfile, chainindex, weights_path='./model.pkl', temperature=298.15, verbose=False
    ):
    """
    Run ProAffinity inference on a given PDB file and chain index.
    
    Args:
        pdbfile (str): Path to the PDB or PDBQT file.
        chainindex (str): List of chains (e.g. "AB,C").
        weights_path (str): Path to the model weights.
        temperature (float): Temperature in Kelvin (default is 298.15) to convert K from.
        
    Returns:
        float: Predicted dG (kJ/mol) value.
    """

    pdbfile = pdbfile
    chainindex = chainindex.split(',')

    # get 2 residue lists from file of one pdb
    # each residue list may contain several res lists, each res list represents a chain

    prot_pair = chainindex

    if len(prot_pair) != 2:
        print('no 2 pro-pro pair!!')

    # Get residues from the PDB file
    if verbose: print(prot_pair)
    residuelistA = get_residue_list_from_file(pdbfile, prot_pair[0])
    residuelistB = get_residue_list_from_file(pdbfile, prot_pair[1])

    adjust_residuelist_num(residuelistA)
    adjust_residuelist_num(residuelistB)

    # Create a sequence for each chain
    if len(residuelistA) == 0 or len(residuelistB) == 0:
        print('no residue list')

    seqA = []

    for reslist in residuelistA: # for each chain in the residuelistA:
        chain_res = []
        for res in reslist:
            chain_res.append(res['type'])
        chain_res = ''.join(chain_res)
        seqA.append(chain_res)

    seqB = []

    for reslist in residuelistB: # for each chain in the residuelistA:
        chain_res = []
        for res in reslist:
            chain_res.append(res['type'])
        chain_res = ''.join(chain_res)
        seqB.append(chain_res)
        
    seqA_len = sum([len(seq) for seq in seqA])

    adjust_residuelist_num(residuelistA)
    adjust_residuelist_num(residuelistB)

    concat_reslistA = list(chain.from_iterable(residuelistA))
    concat_reslistB = list(chain.from_iterable(residuelistB))

    # Generate graphs
    thred = 6
    inter_pairs = get_interaction_residue_pair_new(thred, concat_reslistA, concat_reslistB)
    edge_index = get_edge_index(inter_pairs, seqA_len)
    info_save = (inter_pairs, edge_index, (seqA, seqB))

    # get 2 residue lists from file of one pdb
    # each residue list may contain several res lists, each res list represents a chain

    residuelistA = get_residue_list_from_file(pdbfile, prot_pair[0])
    residuelistB = get_residue_list_from_file(pdbfile, prot_pair[1])

    adjust_residuelist_num(residuelistA)
    adjust_residuelist_num(residuelistB)

    if len(residuelistA) == 0 or len(residuelistB) == 0:
        print('no residue list')

    seqA = []

    for reslist in residuelistA: # for each chain in the residuelistA:
        chain_res = []
        for res in reslist:
            chain_res.append(res['type'])
        chain_res = ''.join(chain_res)
        seqA.append(chain_res)

    seqB = []

    for reslist in residuelistB: # for each chain in the residuelistA:
        chain_res = []
        for res in reslist:
            chain_res.append(res['type'])
        chain_res = ''.join(chain_res)
        seqB.append(chain_res)

    adjust_residuelist_num(residuelistA)
    adjust_residuelist_num(residuelistB)

    concat_reslistA = list(chain.from_iterable(residuelistA))
    concat_reslistB = list(chain.from_iterable(residuelistB))


    thred = 3.5
    intra_pairsA = get_interaction_residue_pair_new_indi(thred, concat_reslistA, concat_reslistA)
    intra_pairsB = get_interaction_residue_pair_new_indi(thred, concat_reslistB, concat_reslistB)

    edge_index1 = get_edge_index_indi(intra_pairsA)
    edge_index2 = get_edge_index_indi(intra_pairsB)

    info_save1 = (intra_pairsA, edge_index1, seqA)
    info_save2 = (intra_pairsB, edge_index2, seqB)
                    
    # Get ESM embedding from the sequences
    # Suppress all warnings during transformers model loading
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from transformers import AutoTokenizer, EsmModel
        import torch

        tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
        model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model.eval()
    torch.set_grad_enabled(False)  

    from torch_geometric.data import Data
    

    atom_type = ['A', 'C', 'OA', 'N', 'NA' 'SA', 'HD']

    atom_pair = ['A_A', 'A_C', 'A_OA', 'A_N', 'A_NA', 'A_SA', 'A_HD', 
                'C_C', 'C_OA', 'C_N', 'C_NA', 'C_SA', 'C_HD',
                'OA_OA', 'OA_N', 'OA_NA', 'OA_SA', 'OA_HD',
                'N_N', 'N_NA', 'N_SA', 'N_HD', 
                'NA_NA', 'NA_SA', 'NA_HD',
                'SA_SA', 'SA_HD',
                'HD_HD']

    bin_number = 10
    type_number = len(atom_pair)
    inter_distance = 6
    intra_distance = 3.5

    seqA, seqB = get_fasta_seq(info_save[2])
    output1_list = []
    output2_list = []

    for fasta in seqA:
        input1 = tokenizer(fasta, return_tensors="pt")
        output1 = model(**input1)
        last_hidden_state1 = output1.last_hidden_state
        last_hidden_state1 = torch.squeeze(last_hidden_state1)
        # get the token from the 2nd to the 2nd last one
        last_hidden_state1 = last_hidden_state1[1:-1]
        output1_list.append(last_hidden_state1)


    for fasta in seqB:
        input2 = tokenizer(fasta, return_tensors="pt")
        output2 = model(**input2)
        last_hidden_state2 = output2.last_hidden_state
        last_hidden_state2 = torch.squeeze(last_hidden_state2)
        # get the token from the 2nd to the 2nd last one
        last_hidden_state2 = last_hidden_state2[1:-1]
        output2_list.append(last_hidden_state2)


    x1 = torch.cat(output1_list, 0)
    x2 = torch.cat(output2_list, 0)
    x = torch.cat((x1, x2), 0)

    pairs = info_save[0]
    edge_index = info_save[1]

    try:
        edge_feature = []
        for pair in pairs:

            edge_encoding = np.zeros(type_number * bin_number)
            residueA = pair[0]
            residueB = pair[1]

            for atom1 in residueA['atoms']:
                x1 = float(atom1['x'])
                y1 = float(atom1['y'])
                z1 = float(atom1['z'])
                type1 = atom1['pdbqt_type']

                for atom2 in residueB['atoms']:
                    x2 = float(atom2['x'])
                    y2 = float(atom2['y'])
                    z2 = float(atom2['z'])
                    type2 = atom2['pdbqt_type']

                    dis = get_distance(x1, y1, z1, x2, y2, z2)

                    bin_n = math.ceil(dis / (inter_distance / bin_number))
                    #print(bin_n)
                    if bin_n > 10:
                        bin_n = 10

                    if type1 + '_' + type2 in atom_pair:
                        pair_type = type1 + '_' + type2       
                    elif type2 + '_' + type1 in atom_pair:
                        pair_type = type2 + '_' + type1
                    else:
                        print('no match type!' + 'type1:' + type1 + 'type2:' + type2)
                        pair_type = 'others'
                        
                    pair_type_index = atom_pair.index(pair_type)
                    encoding_index = (bin_n - 1) * type_number + pair_type_index
                    edge_encoding[encoding_index] = edge_encoding[encoding_index] + 1

            edge_feature.append(torch.from_numpy(edge_encoding))                

        edge_feature = torch.stack(edge_feature, 0) 
        edge_feature = torch.cat((edge_feature, edge_feature), 0)

    except Exception as e:
        print(e)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_feature)

    # for individual graph

    seq1 = info_save1[2]

    output1_list = []

    for fasta in seq1:
        input1 = tokenizer(fasta, return_tensors="pt")
        output1 = model(**input1)
        last_hidden_state1 = output1.last_hidden_state
        last_hidden_state1 = torch.squeeze(last_hidden_state1)
        # get the token from the 2nd to the 2nd last one
        last_hidden_state1 = last_hidden_state1[1:-1]
        output1_list.append(last_hidden_state1)

    x_indi_1= torch.cat(output1_list, 0)

    pairs = info_save1[0]
    edge_index = info_save1[1]

    try:
        edge_feature = []
        for pair in pairs:

            edge_encoding = np.zeros(type_number * bin_number)
            residueA = pair[0]
            residueB = pair[1]

            for atom1 in residueA['atoms']:
                x1 = float(atom1['x'])
                y1 = float(atom1['y'])
                z1 = float(atom1['z'])
                type1 = atom1['pdbqt_type']

                for atom2 in residueB['atoms']:
                    x2 = float(atom2['x'])
                    y2 = float(atom2['y'])
                    z2 = float(atom2['z'])
                    type2 = atom2['pdbqt_type']

                    dis = get_distance(x1, y1, z1, x2, y2, z2)

                    bin_n = math.ceil(dis / (intra_distance / bin_number))
                    #print(bin_n)
                    if bin_n > 10:
                        bin_n = 10

                    if type1 + '_' + type2 in atom_pair:
                        pair_type = type1 + '_' + type2       
                    elif type2 + '_' + type1 in atom_pair:
                        pair_type = type2 + '_' + type1
                    else:
                        print('no match type!' + 'type1:' + type1 + 'type2:' + type2)
                        pair_type = 'others'
                        
                    pair_type_index = atom_pair.index(pair_type)
                    encoding_index = (bin_n - 1) * type_number + pair_type_index
                    edge_encoding[encoding_index] = edge_encoding[encoding_index] + 1

            edge_feature.append(torch.from_numpy(edge_encoding))                

        edge_feature = torch.stack(edge_feature, 0) 

    except Exception as e:
        print(e)

    data1 = Data(x=x_indi_1, edge_index=edge_index, edge_attr=edge_feature)

    seq2 = info_save2[2]

    output2_list = []

    for fasta in seq2:
        input2 = tokenizer(fasta, return_tensors="pt")
        output1 = model(**input2)
        last_hidden_state2 = output2.last_hidden_state
        last_hidden_state2 = torch.squeeze(last_hidden_state2)
        # get the token from the 2nd to the 2nd last one
        last_hidden_state2 = last_hidden_state2[1:-1]
        output2_list.append(last_hidden_state2)

    x_indi_2= torch.cat(output2_list, 0)

    pairs = info_save2[0]
    edge_index = info_save2[1]

    try:
        edge_feature = []
        for pair in pairs:

            edge_encoding = np.zeros(type_number * bin_number)
            residueA = pair[0]
            residueB = pair[1]

            for atom1 in residueA['atoms']:
                x1 = float(atom1['x'])
                y1 = float(atom1['y'])
                z1 = float(atom1['z'])
                type1 = atom1['pdbqt_type']

                for atom2 in residueB['atoms']:
                    x2 = float(atom2['x'])
                    y2 = float(atom2['y'])
                    z2 = float(atom2['z'])
                    type2 = atom2['pdbqt_type']

                    dis = get_distance(x1, y1, z1, x2, y2, z2)

                    bin_n = math.ceil(dis / (intra_distance / bin_number))
                    #print(bin_n)
                    if bin_n > 10:
                        bin_n = 10

                    if type1 + '_' + type2 in atom_pair:
                        pair_type = type1 + '_' + type2       
                    elif type2 + '_' + type1 in atom_pair:
                        pair_type = type2 + '_' + type1
                    else:
                        print('no match type!' + 'type1:' + type1 + 'type2:' + type2)
                        pair_type = 'others'
                        
                    pair_type_index = atom_pair.index(pair_type)
                    encoding_index = (bin_n - 1) * type_number + pair_type_index
                    edge_encoding[encoding_index] = edge_encoding[encoding_index] + 1

            edge_feature.append(torch.from_numpy(edge_encoding))                

        edge_feature = torch.stack(edge_feature, 0) 

    except Exception as e:
        print(e)


    data2 = Data(x=x_indi_2, edge_index=edge_index, edge_attr=edge_feature)

    graph = data
    graph1 = data1
    graph2 = data2
    graph.edge_attr = graph.edge_attr.float()

    graph1.edge_attr = [attr.float() for attr in graph1.edge_attr]
    graph2.edge_attr = graph2.edge_attr.float()

    datalist_inter = []
    datalist_intra1 = []
    datalist_intra2 = []
    datalist_inter.append(graph)
    datalist_intra1.append(graph1)
    datalist_intra2.append(graph2)

    test_loader_inter = DataLoader(datalist_inter, batch_size=1)
    test_loader_intra1 = DataLoader(datalist_intra1, batch_size=1)
    test_loader_intra2 = DataLoader(datalist_intra2, batch_size=1)

    in_channels = data.num_node_features
    hidden_channels = 256
    out_channels = 64
    linear_out1 = 32   
    linear_out2 = 1
    edge_dim = data.num_edge_features
    num_layers = 3 
    num_timesteps = 2
    dropout = 0.5

    devicename = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(devicename)
    model = GraphNetwork(in_channels, hidden_channels, out_channels, edge_dim, num_layers, num_timesteps, dropout, linear_out1, linear_out2).to(device)

    state_dict = torch.load(weights_path, map_location=devicename)

    new_state_dict = {}
    for key, value in state_dict.items():
        new_key = key.replace("lin_src", "lin").replace("lin_dst", "lin")
        new_state_dict[new_key] = value

    model.load_state_dict(new_state_dict, strict=False)


    # Assuming model is your GNN model and dataloader is your test dataloader
    model.eval()  # Set the model to evaluation mode
    

    all_predictions = []

    with torch.no_grad():
        # Disable gradient computation during testing
        for batch_inter, batch_intra1, batch_intra2 in zip(test_loader_inter, test_loader_intra1, test_loader_intra2):
            # Assuming batch contains input data 'x' and true values 'y_true'
            batch_inter = batch_inter.to(device)
            batch_intra1 = batch_intra1.to(device)
            batch_intra2 = batch_intra2.to(device)
            if isinstance(batch_intra1.edge_attr, list):
                batch_intra1.edge_attr = torch.stack(batch_intra1.edge_attr, dim=0).to(device) 
            if isinstance(batch_intra2.edge_attr, list):
                batch_intra2.edge_attr = torch.stack(batch_intra2.edge_attr, dim=0).to(device)        # y_true = batch_inter.y
            
            # Get model predictions for the current batch
            y_pred = model(batch_inter, batch_intra1, batch_intra2)
            y_pred = torch.squeeze(y_pred)
            
            # Store predictions and true values
            all_predictions.append(y_pred.cpu().numpy())

    for i in range(len(all_predictions)):
        # Check if the current item is a scalar by examining its dimensionality
        if all_predictions[i].ndim == 0:
            # Convert scalar to a 1D array and update the item in the list
            all_predictions[i] = np.array([all_predictions[i]])

    # Concatenate all predictions and true values
    all_predictions = np.concatenate(all_predictions, axis=0)
    K = 10**(all_predictions[0])
    R = 8.314 / 1000 # kJ/(mol*K)
    dG = -R * temperature * np.log(K)
    return dG

# Example usage
if __name__ == "__main__": 
    run_proaffinity_inference("./pdbfiles/1i4d.pdbqt", "D,AB", verbose=True)


