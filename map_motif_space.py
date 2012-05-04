#
# Input: a directory of mlib files
#
# Output: a graph where node size = # of motifs bounds by a PWM, edge weight = # of motifs jointly bound by two PWMs
#

from argparser import *
import matplotlib.pyplot as plt
import networkx as nx
import os

ap = ArgParser(sys.argv)

def build_mlib_hash(genes_files):
    """genes_files[gene] = path to mlib"""
    """Returns ret[gene] = list of motifs"""
    ret = {}
    for gene in genes_files.keys():
        ret[gene] = []
        f = genes_files[gene]
        fin = open(f, "r")
        lines = fin.readlines()
        fin.close()
        for l in lines:
            if l.__len__() > 2 and False == l.startswith("#"):
                ret[gene].append( l.strip() )
        print gene, ret[gene]
    return ret
                

def get_mlib_files(dirpath):
    """Input: directory path, output = list of mlib files."""
    mlib_files = {}
    for f in os.listdir( dirpath ):
        if f.__contains__("mlib"):
            tokens = f.split(".")
            gene = tokens[1]
            mlib_files[gene] = dirpath + "/" + f
    return mlib_files

def plot_mlib_distribution(tf_m):
    mliblens = []
    for tf in tf_m.keys():
        mliblens.append( tf_m[tf].__len__() )
    plt.hist(mliblens, 20)
    plt.show()

def intersect(a, b):
     return list(set(a) & set(b))

def plot_motif_space(tf_m):
    G = nx.Graph()
    for tf in tf_m.keys():
        G.add_node(tf, size=1.0*tf_m[tf].__len__())
    
    tfs = tf_m.keys()
    for i in range(0, tfs.__len__()):
        for j in range(i+1, tfs.__len__()):
            x = intersect(tf_m[ tfs[i] ], tf_m[ tfs[j] ]).__len__()
            if x > 0:
                print tfs[i], tfs[j], x
                G.add_edge(tfs[i], tfs[j], weight=1.0*x)
            
    plt.figure(figsize=(8,8))
    pos=nx.spring_layout(G,iterations=100)
    nodesize=[]
    for v in G.node:
        nodesize.append(G.node[v]["size"])
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color="blue", alpha=0.5, linewidths=0.1)

    for e in G.edges():
        print e
        edgewidth = [ G.get_edge_data(e[0],e[1])["weight"] ]
        this_edge = [ e ]
        print this_edge, edgewidth
        #print [(pos[e[0]],pos[e[1]]) for e in this_edge]
        nx.draw_networkx_edges(G, pos, edgelist = this_edge, width = edgewidth)

    nx.draw_networkx_labels(G, pos, font_size=9, font_family="Helvetica")
    
    plt.show()


#
#
# MAIN:
#
#
mlib_dir = ap.getArg("--mlibdir")
output_dir = ap.getArg("--outputdir")
mlib_files = get_mlib_files(mlib_dir)
tf_m = build_mlib_hash(mlib_files)
#plot_mlib_distribution( tf_m )
plot_motif_space( tf_m )