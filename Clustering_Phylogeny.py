from Bio.Cluster import kcluster 
from Bio.Cluster import kmedoids
from Bio.Cluster import distancematrix
from Bio.SeqIO.XdnaIO import _parse_feature_description
import numpy as np
from Bio import SeqIO
from Bio import Phylo
import matplotlib as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os


def read_seq(seq_file_path, file_format):
    '''A method to read the sequences from file'''

    #Read MSA and convert it to a dictionary
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_file_path, file_format)) 

    #Save Sequences and sequence ID
    sequences = []

    with open(seq_file_path) as handle:
        for record in SeqIO.parse(handle, file_format):
            sequences.append(str(record.seq))
    return sequences, record_dict

def cluster_sequences(seq_file_path, file_format):
    '''A method for clustering sequences'''
    #Clustering my Sequences
    sequences = read_seq(seq_file_path, file_format)[0]
    sequence_dict = read_seq(seq_file_path, file_format)[1]
    seq_matrix = np.asarray([np.fromstring(s, dtype=np.uint8) for s in sequences])
    # clusterid, error,found = kcluster(seq_matrix, nclusters=8, npass=4) 
    # print(clusterid) 
    # print(found) 
    # print(error)
    matrix = distancematrix(seq_matrix, mask=None, weight=None, transpose=0, dist='e')
    medroid_clusterid, error, nfound = kmedoids(matrix, nclusters=8, npass=4,initialid=None)
    print("Kmedroid")
    print(medroid_clusterid) 
    print(nfound)
    zipped_seq = zip(medroid_clusterid,sequences)
    cluster_dict = {}
    for (key, value) in zipped_seq:
        for id, value1  in sequence_dict.items():
            if value == value1.seq:
                description = value1.description.split(" ",2)
                value1.id =description[2]
                value1.name = description[1]
                cluster_dict.setdefault(key, []).append(value1)
    for key,value in cluster_dict.items():
        filename= "aligned_sequences/Cluster_files"+str(key) + ".fa"
        SeqIO.write(value, filename, "fasta")
    return cluster_dict


def build_cluster_tree(cluster_dict):
    '''Building cluster  sequence tree'''
    for key,value in cluster_dict.items():
        align = MultipleSeqAlignment(value)
        #Protein calculator with ‘blosum62’ model:
        calculator = DistanceCalculator('blosum62')
        dm = calculator.get_distance(align)
        #Constructor to construct tree
        constructor = DistanceTreeConstructor()
        #Constructs and returns an Unweighted Pair Group Method with Arithmetic mean (UPGMA) tree.
        upgmatree = constructor.upgma(dm)
        #print(upgmatree)

        #Make tree Parsimonius
        scorer = Phylo.TreeConstruction.ParsimonyScorer()
        searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
        constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, upgmatree)
        pars_tree = constructor.build_tree(align)
        #print(pars_tree)
        
        #write tree to file
        treename = "Trees/Cluster/"+str(key)+".xml"
        Phylo.write(pars_tree, treename, "phyloxml")
        draw_tree(pars_tree)
        
    

def build_tree(aligned_sequence, aligned_squence_filetype, seq_type):
    '''A methid for building trees using aligned sequences'''
    #loading alignment
    aln = AlignIO.read(open(aligned_sequence), aligned_squence_filetype)
    #print(aln)
    for record in aln:
        description = record.description.split(" ",2)
        record.id = description[2]
        record.name = description[1]
        

    if seq_type == "dna":
        #DNA calculator with ‘identity’ model:
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        #print(dm)
    elif seq_type == "protein":
        #Protein calculator with ‘blosum62’ model:
        calculator = DistanceCalculator('blosum62')
        dm = calculator.get_distance(aln)
        #print(dm)
    else:
        print("please pass in right sequence type e.g. protein or dna")

    #Constructor to construct tree
    constructor = DistanceTreeConstructor()
    #Constructs and returns an Unweighted Pair Group Method with Arithmetic mean (UPGMA) tree.
    upgmatree = constructor.upgma(dm)
    #print(upgmatree)

    #Make tree Parsimonius
    scorer = Phylo.TreeConstruction.ParsimonyScorer()
    searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
    constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, upgmatree)
    pars_tree = constructor.build_tree(aln)
    print(pars_tree)
    # ##Construct and return a Neighbor Joining tree.
    # # njtree = constructor.nj(dm)
    # # print(njtree)
    return pars_tree






def consensus_tree(tree_file_path):
    '''This method is to create Consensus tree'''
    tree_list = []
    for fname in os.listdir(tree_file_path):
        if fname.startswith != tree_file_path+"all-alignment":
            tree = Phylo.read(fname, "phyloxml")
            tree_list.append(tree)  
    strict_consensus = Phylo.Consensus.strict_consensus(tree_list)
    print(strict_consensus)
    majority_consensus = Phylo.Consensus.majority_consensus(tree_list, cutoff=0.5)
    print(majority_consensus)
    adam_consensus = Phylo.Consensus.adam_consensus(tree_list)
    print(adam_consensus)
    return strict_consensus, majority_consensus, adam_consensus

def draw_tree(tree):
    '''A method for visualizing trees through drawings'''
    #draw tree
    tree = tree.as_phyloxml()
    tree.root.color = (128, 128, 128)
    Phylo.draw(tree) 

    ##Phylogenetics - visualizing tree
    # tree = Phylo.read("simple.dnd", "newick")
    # print(tree)
    # #Phylo.draw_ascii(tree)
    # tree.rooted = True
    # #Phylo.draw(tree)
    # tree = tree.as_phyloxml()
    # tree.root.color = (128, 128, 128)
    # mrca = tree.common_ancestor({"name": "E"}, {"name": "F"})
    # mrca.color = "salmon"
    # tree.clade[1, 0].color = "blue"
    # #Phylo.draw(tree)



def main_program():
    '''Method for running the full program'''
    for fname in os.listdir('aligned_sequences'):
        if fname.endswith('.fa'):
            filepath = 'aligned_sequences/'+fname
            seq_file_type = "fasta"
            built_tree = build_tree(filepath,seq_file_type,"protein")
            #write tree to file
            treename = "Trees/"+fname.split(".")[0]+".xml"
            Phylo.write(built_tree, treename, "phyloxml")
            draw_tree(built_tree)
            


cluster_seq = cluster_sequences("aligned_sequences/all-alignment.fa","fasta")
build_cluster_tree(cluster_seq)
main_program()
consensus_tree = consensus_tree("Trees")
cluster_consensus_tree = consensus_tree("Trees/Cluster")
draw_tree(consensus_tree[0])
draw_tree(consensus_tree[1])
draw_tree(consensus_tree[2])
draw_tree(cluster_consensus_tree[0])
draw_tree(cluster_consensus_tree[1])
draw_tree(cluster_consensus_tree[2])