from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
import numpy as np

# TODO change atom list to all atoms
_possible_atom_list =  ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg',
    'Na', 'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn',
    'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 
    'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt','Hg', 'Pb', 'As', 'UNK']
_possible_numH_list = [0, 1, 2, 3, 4]
_possible_valence_list = [0, 1, 2, 3, 4, 5, 6]
_possible_formal_charge_list = [-3, -2, -1, 0, 1, 2, 3]
_possible_degree_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
_possible_hybridization_list = [
    Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
    Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
    Chem.rdchem.HybridizationType.SP3D2]
_possible_number_radical_e_list = [0, 1, 2]
_possible_chirality_list = ['R', 'S']


def one_hot_encoding(x, set):
    if x not in set:
        raise Exception("input {0} not in allowable set{1}:".format(x, set))
    return list(map(lambda s: x == s, set))


def one_hot_encoding_unk(x, set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in set:
        x = set[-1]
    return list(map(lambda s: x == s, set))


def encode_atom(atom, bool_id_feat=False, 
                explicit_H=False, use_chirality=False):
    """
    From deepchem.feat.graph_features
    """

    # why not one-hot get_formal_charge and get_num_radical_electrons?      
    result = \
        one_hot_encoding_unk(
            atom.GetSymbol(), _possible_atom_list) + \
        one_hot_encoding(atom.GetDegree(), _possible_degree_list) + \
        one_hot_encoding_unk(
            atom.GetImplicitValence(), _possible_valence_list) + \
        [atom.GetFormalCharge(), atom.GetNumRadicalElectrons()] + \
        one_hot_encoding_unk(
            atom.GetHybridization(), _possible_hybridization_list) + \
        [atom.GetIsAromatic()]
            
    if not explicit_H:
        result = result + one_hot_encoding_unk(
            atom.GetTotalNumHs(), _possible_numH_list)
    if use_chirality:
        try:
            result = result + one_hot_encoding_unk(
                atom.GetProp('_CIPCode'), _possible_chirality_list) + \
                [atom.HasProp('_Ch_possible_numH_list = iralityPossible')]
        except:
            result = result + [False, False] + \
                [atom.HasProp('_ChiralityPossible')]
                
    return np.array(result)


def build_graph_from_molecule(mol, use_master_atom=False):
    """
    From deepchem.ConvMolFeaturizer._featurize
    Input:
        rdkit.Chem.rdchem.Mol
        
    Output:
        nodes - np.ndarray of shape (num_atoms, num_feat)
        canon_adj_list - list. index corresponds to the index of node 
                         and canon_adj_list[index] corresponds to indices
                         of the nodes that node i is connected to. 
    """
    idx_nodes = [(atom.GetIdx(), encode_atom(atom)) 
                 for atom in mol.GetAtoms()]
    idx_nodes.sort()
    _, nodes = list(zip(*idx_nodes))
    
    nodes = np.vstack(nodes)

    # Master atom is the "average" of all atoms that is connected to all atom
    # Introduced in https://arxiv.org/pdf/1704.01212.pdf
    if use_master_atom: 
        master_atom_features = np.expand_dims(np.mean(nodes, axis=0), axis=0)
        nodes = np.concatenate([nodes, master_atom_features], axis=0)
        
    edge_list = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) 
                for bond in mol.GetBonds()]
    
    canon_adj_list = [[] for _ in range(len(nodes))] 
    
    for edge in edge_list:
        canon_adj_list[edge[0]].append(edge[1])
        canon_adj_list[edge[1]].append(edge[0])
        
    if use_master_atom:
        fake_atom_index = len(nodes) - 1
        
        for i in range(len(nodes) - 1):
            canon_adj_list[i].append(fake_atom_index)
        
    return (nodes, canon_adj_list)


def smiles2graph(smiles_arr):
    """
    from deepchem.data.data_loader.featurize_smiles_df
    """
    features = []
    invalid_ind = []

    for ind, smiles in enumerate(smiles_arr):
        try:
            mol = Chem.MolFromSmiles(smiles)

            # what are the two lines below doing?
            # Answer found in deepchem.data.data_loader featurize_smiles_df
            # TODO (ytz) this is a bandage solution to reorder the atoms so
            # that they're always in the same canonical order. Presumably this
            # should be correctly implemented in the future for graph mols.
            new_order = rdmolfiles.CanonicalRankAtoms(mol)
            mol = rdmolops.RenumberAtoms(mol, new_order)

            feature = build_graph_from_molecule(mol)
            features.append(feature)

        except:
            invalid_ind.append(ind)

    return features, invalid_ind
