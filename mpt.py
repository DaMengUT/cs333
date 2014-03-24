#!/usr/bin/env python

import random
import newick
from math import log
from pygraph.algorithms.accessibility import connected_components
from pygraph.classes.exceptions import AdditionError
from union_find import *
from util import *

#import sys
#sys.path.append('..')
#sys.path.append('/usr/lib/graphviz/python/')
#sys.path.append('/usr/lib64/graphviz/python/')

def heuristic(sequences, weights):
    def gen_typical_sequence(sequences, seq_names):
        return ''.join([random.choice(x) for x in zip(*[sequences[n] for n in seq_names])])
    result = graph()
    result.add_nodes(sequences.keys())
    internal_node_label_phases = "abcdefghijklmnopqrstuvwxyz"
    root = None
    for phase in internal_node_label_phases:
        gr = heuristic_1(sequences, weights)
        print_dot(result,phase + ".png")
        cc = connected_components(gr)
        #print cc
        size_cc = len(cc)
        if size_cc == 1:
            break

        cc_seqs = {}
        # Transform from {seq_name : component_index}
        # to {component_index : [seq_name]} joining by
        # component_index
        for seq_name in cc.keys():
            if cc[seq_name] not in cc_seqs:
                cc_seqs[cc[seq_name]] = []
            cc_seqs[cc[seq_name]].append(seq_name)

        typical_sequences = {}
        for k in cc_seqs.keys():
            typical_sequences[phase + str(k)] = gen_typical_sequence(sequences, cc_seqs[k])

        # add edges to the result graph, giving each component its own subtree
        for k in cc_seqs.keys():
            result.add_node(phase+str(k))
            for seq in cc_seqs[k]:
                result.add_edge((phase + str(k),seq))

        sequences = typical_sequences
        weights = gen_sequence_weights(sequences)
        root = phase + str(1)
    return (result,root)


def heuristic_1(sequences, weights):
    # find matchings
    gr = graph()
    gr.add_nodes(sequences.keys())
    mins = {}

    def cost(s1,s2):
        if weights.has_key((s1,s2)):
            return weights[(s1,s2)]
        else:
            return float('+inf')
    for s in sequences.keys():
        mins[s] = min([(cost( s,t ),t) for t in sequences.keys()])
        #print mins[s]
        try:
            gr.add_edge((s,mins[s][1]))
        except AdditionError:
            pass
    #print mins
    return gr

def heuristic_3(vertex_set,weights):

    tree = []
    pairs = weights.items()
    uf = UnionFind()
    uf.insert_objects(vertex_set)
    sorted_pairs = sorted(pairs, key=lambda pair: pair[1])
    edge,weight = sorted_pairs[0]
    tree.append(edge[0])
    tree.append(edge[1])
    sorted_pairs.remove((edge,weight))

    #print vertex_set
    tree = get_sequence(tree, sorted_pairs)
    #print tree
    return tree

def get_sequence(tree, sorted_pairs):
    last_sequence_added = tree[-1]
    min_weight = 1000
    save = ("","")
    for edge,weight in sorted_pairs:
        if edge[0] == last_sequence_added and not edge[1] in tree and (weight < min_weight):
            save = edge
            min_weight = weight
        if edge[1] == last_sequence_added and not edge[0] in tree and (weight < min_weight):
            save = edge
            min_weight = weight

    if save == ("",""):
        return tree

    if save[0] == last_sequence_added:
        tree.append(save[1])
    else:
        tree.append(save[0])

    if save != ("",""):
        return get_sequence(tree, sorted_pairs)
    else:
        return tree


def heuristic_2(tree):
    gg = graph()
    gg.add_nodes(tree)

    tree_depth = int(log(len(tree), 2)) + 1
    seqs=tree
    for _,phase in zip(range(tree_depth),"abcdefghijklmnopqrstuvwxyz"):
        new_seqs=[]
        for i,(x,y) in enumerate(grouped(seqs, 2)):
            gg.add_node(phase + str(i))
            new_seqs.append(phase + str(i))
            gg.add_edge((phase + str(i), x))
            gg.add_edge((phase + str(i), y))
        if len(new_seqs) == 0:
            break
        if len(seqs) % 2 == 1:
            # There's a sequence left over
            gg.add_edge((new_seqs[-1], seqs[-1]))
        seqs=new_seqs

    print_dot(gg, "asd" + ".png")
    return (gg, tree[0])

def gen_sequence_weights (strings):
    """ Takes a dictionary of sequence names and sequences """
    distances = {}
    for x,y in product(strings.keys(), strings.keys()):
        if (x!=y):
            distances[(x,y)] = hamming_distance(strings[x], strings[y])
    return distances

def prune():
    pass

data = []

warnow_files = ["true.100.fasta", "true.25.fasta"]
phyML_files = ["all_seq"]


for k in range(1, 5000):
    strings = read_phyML(phyML_files[0],k)
    distances = gen_sequence_weights(strings)
    true_tree = read_phyML_tree("all_tree",k)
    tt_cost = calc_fitch_cost(true_tree,strings,1)

    (score,tree) = kruskals(strings.keys(),distances)

    h1,start = heuristic(strings, distances)
    h1_cost = calc_fitch_cost(h1,strings,start)

    h3_tree = heuristic_3(strings.keys(),distances)
    h2,start2 = heuristic_2(h3_tree)
    h2_cost = calc_fitch_cost(h2,strings,start2)


    data.append([float(c) / tt_cost for c in (tt_cost,score,h1_cost,h2_cost)])
    s=data[k-1]
    print ",".join([str(x) for x in s])
