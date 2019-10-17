import networkx as nx
from rdkit.Chem.rdmolfiles import MolFromPDBFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def mol_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic())

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())

    return G

protein = MolFromPDBFile('./data/pdbbind/v2018/10gs/10gs_protein_fixed.pdb')
protein_graph = mol_to_nx(protein)
pos = nx.spring_layout(protein_graph)

f = plt.figure()
nx.draw(protein_graph, pos, width=1, node_size=1)
f.savefig('graph.png')
