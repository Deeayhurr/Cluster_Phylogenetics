from Bio.Cluster import kcluster 
from Bio.Cluster import kmedoids
from Bio.Cluster import distancematrix
from Bio.SeqIO.XdnaIO import _parse_feature_description
import numpy as np
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo import Consensus
import matplotlib as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.PhyloXML import Phylogeny
from Bio.Align import MultipleSeqAlignment
import os, shutil
import glob



def read_seq(seq_file_path, file_format):
    '''A method to read the sequences from file'''

    print("reading Msa file sequence")

    #Read MSA and convert it to a dictionary
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_file_path, file_format)) 

    #Save Sequences and sequence ID
    sequences = []

    with open(seq_file_path) as handle:
        for record in SeqIO.parse(handle, file_format):
            sequences.append(str(record.seq))
    print("reading Msa file sequence completed")
    return sequences, record_dict

def cluster_sequences(seq_file_path, file_format):
    '''A method for clustering sequences'''

    print("starting: clustering  sequence")
    #Clustering my Sequences
    sequences = read_seq(seq_file_path, file_format)[0]
    sequence_dict = read_seq(seq_file_path, file_format)[1]
    seq_matrix = np.asarray([np.fromstring(s, dtype=np.uint8) for s in sequences])
    # clusterid, error,found = kcluster(seq_matrix, nclusters=8, npass=4) 
    # print(clusterid) 
    # print(found) 
    # print(error)
    matrix = distancematrix(seq_matrix, mask=None, weight=None, transpose=0, dist='e')
    medroid_clusterid, error, nfound = kmedoids(matrix, nclusters=8, npass=3,initialid=None)
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
        filename= "aligned_sequences/Cluster_files/"+str(key) + ".fa"
        SeqIO.write(value, filename, "fasta")
    print("ended: clustering  sequence")
    return cluster_dict


def build_cluster_tree(cluster_dict):
    '''Building cluster  sequence tree'''
    print("starting:  building clustering trees")
    cluster_tree_lists = []
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
        cluster_tree_lists.append(pars_tree)
        
        #write tree to file
        treename = "Trees/Cluster/"+str(key)+".xml"
        Phylo.write(pars_tree, treename, "phyloxml")
        #draw_tree(pars_tree)
    print("ended:  building clustering trees")
    return cluster_tree_lists
 

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
    if not aligned_sequence.endswith("all-alignment.fa"):
        #Make tree Parsimonius
        scorer = Phylo.TreeConstruction.ParsimonyScorer()
        searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
        constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, upgmatree)
        pars_tree = constructor.build_tree(aln)
        #print(pars_tree)

        return pars_tree
    else:
        return upgmatree
    # ##Construct and return a Neighbor Joining tree.
    # # njtree = constructor.nj(dm)
    # # print(njtree)
    
def group_consensus_tree(trees):
    print("starting:  building group consensus trees")
    majority_consensus = Consensus.majority_consensus(trees, cutoff=0.5)
    #strict__consensus = Consensus.strict_consensus(trees)
    # adam__consensus = Consensus.adam_consensus(trees)
    print("ended:  building group consensus trees")
    return majority_consensus

def consensus_tree(aligned_sequence, aligned_squence_filetype):
    '''This method is to create Consensus tree'''
    aln = None
    try : 
        aln = AlignIO.read(open(aligned_sequence), aligned_squence_filetype)
    except Exception as e:
        print (e.__cause__, aligned_sequence)
        exit()
    for record in aln:
        description = record.description.split(" ",2)
        record.id = description[2]
        record.name = description[1]

    calculator = DistanceCalculator("blosum62")
    constructor = DistanceTreeConstructor(calculator)
    consensus = Consensus.bootstrap_consensus(aln, 100, constructor, Consensus.majority_consensus)
    return consensus


def draw_tree(tree):
    '''A method for visualizing trees through drawings'''
    #draw tree
    #tree = Phylo.read("Trees/all-alignment.xml", "phyloxml")
    #tree = Phylogeny.from_tree(tree)
    tree = tree.as_phyloxml()
    tree.root.color = "black"

    all_protein_names = {"albumin" : "red", "alpha" : "blue", "globin" : "yellow", "glucagon" : "orange", "glyco" : "green", "myoglobin" : "purple", "Sirtuin" : "pink", "ZP2" : "cyan", "zona" : "cyan"}
    all_animal_names = {"Papio anubis" : "brown", "Hyaena hyaena" : "gold", "Rhinopithecus roxellana" : "fuchsia", "Piliocolobus tephrosceles" : "lime", "Pongo abelii" : "magenta", "Gorilla gorilla gorilla" : "maroon", "Pan paniscus" : "aqua", "Homo sapiens" : "navy"}

    
    for key, value in all_protein_names.items():
        name = []
        for clade in tree.find_clades({"name" : f".*{key}.*"}):
            name.append(clade.name)
        protein_name = tree.common_ancestor(name)
        protein_name.color = value
                
     
    for key, value in all_animal_names.items():
        for clade in tree.find_clades({"name" : f".*{key}.*"}):
            animal_name=clade.name
            animal_name = tree.common_ancestor(animal_name)
            animal_name.color = value


    #animal_match = tree.find_elements(name='hyaena hyaena')
    #animal_match.color = "red"
    Phylo.draw(tree)  
    ##Phylogenetics - visualizing tree
    # tree = Phylo.read("simple.dnd", "newick")
    # print(tree)
    # #Phylo.draw_ascii(tree)
    # tree.rooted = True
    # tree.clade[1, 0].color = "blue"
    # #Phylo.draw(tree)


def main_program():
    '''Method for running the full program'''
    #build cluster tree
    cluster_seq = cluster_sequences("aligned_sequences/all-alignment.fa","fasta")
    cluster_tree_lists = build_cluster_tree(cluster_seq)

    folder = "aligned_sequences/Cluster_files"
    for file in glob.glob(f'{folder}/*'):
        os.remove(file)
    folder = "Trees/Cluster"
    for file in glob.glob(f'{folder}/*'):
        if os.path.isfile(file):
            os.remove(file)

    #build normal and consensus (cluster and protein)tree for all groups
    print("starting:  building protein and consensus trees")
    trees_list = []
    consensus_tree_list = []
    for fname in os.listdir('aligned_sequences'):
        if fname.endswith('.fa'):
            filepath = 'aligned_sequences/'+fname
            seq_file_type = "fasta"
            built_tree = build_tree(filepath,seq_file_type,"protein")
            #write tree to file
            treename = "Trees/"+fname.split(".")[0]+".xml"
            Phylo.write(built_tree, treename, "phyloxml")
            if not fname.startswith("all-alignment"):
                trees_list.append(built_tree)
                consens_tree = consensus_tree(filepath,seq_file_type)
                consensus_treename = "Trees/Consensus_trees/"+fname.split(".")[0]+".xml"
                Phylo.write(consens_tree, consensus_treename, "phyloxml")
            else:
                all_tree = built_tree
        elif fname == "Cluster_files":
            for filename in os.listdir('aligned_sequences/'+fname):
                if filename.endswith('.fa'): 
                    seq_file_type = "fasta"
                    consens_tree = consensus_tree('aligned_sequences/'+fname+"/"+filename,seq_file_type)
                    consensus_treename = "Trees/Consensus_trees/Cluster/"+filename.split(".")[0]+".xml"
                    Phylo.write(consens_tree, consensus_treename, "phyloxml")   
                #draw_tree(built_tree)
        consensus_tree_list.append(consens_tree)
    print("ended:  building protein and consensus trees")

    print(f'trees_list is {trees_list} having length{len(trees_list)}')
    print(f'cluster trees_list is {cluster_tree_lists} having length{len(cluster_tree_lists)}')
    protein_tree= group_consensus_tree(trees_list)
    cluster_tree = group_consensus_tree(cluster_tree_lists)
    
    #draw group consensus tree protein
    print("starting:  drawing  trees")
    print("drawing group consensus tree protein")
    draw_tree(protein_tree)
    
    #draw group consensus tree cluster
    print("drawing group consensus tree cluster")
    draw_tree(cluster_tree)

    #draw consensus tree for 1 cluster, 1protein and all protein
    draw_tree(trees_list[1])
    draw_tree(consensus_tree_list[2])
    draw_tree(consensus_tree_list[3])
    draw_tree(all_tree)
    print("ended: drawing trees")

    



main_program()
 
# tree=""
# draw_tree(tree)
