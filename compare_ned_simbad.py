import numpy as np
from treelib import Tree, Node
from astropy.table import Table
import string

def do_compare(intab, simbad_file):
    # read data - a table with objects appearing in both NED and SIMBAD
    data_tab = Table.read(intab)
    # the 2 things we care about are the NED and SIMBAD classes
    nclass = data_tab['']
    sclass = data_tab['']
    
    # read SIMBAD hierarchy
    simbad_tree = tree_read(simbad_file)

    # compute match_score
    match_score = match_flat_tree(nclass,sclass,simbad_tree)

    # do something interesting with the output
    # like plot it
    # or save it

    return

EXACT_MATCH = 0
ONE_UP = 1
ONE_DOWN = 2
NON_MATCH = 5
NOT_IN_TREE = 10

def match_flat_tree(flatclass, treeclass, tree, tree_dict):
    '''match entries between NED and SIMBAD classes'''
    if len(flatclass) != len(treeclass):
        print('flatclass and treeclass must have same length!')
        return

    score = np.empty(len(flatclass))
    for i in range(0, len(flatclass)):
#        print('Tree class is leaf:',tree[tree_dict[treeclass[i]]].is_leaf())
#        print('Tree class is root:',tree_dict[treeclass[i]] == tree.root)
        if flatclass[i] not in tree_dict: # flat identifier is not in tree, as node or leaf
            score[i] = NOT_IN_TREE
        elif flatclass[i] == treeclass[i]: # exact match
            score[i] = EXACT_MATCH
        elif tree_dict[treeclass[i]] != tree.root:
            if flatclass[i] == tree.parent(tree_dict[treeclass[i]]).tag:  # flat class is one level up from tree class
                score[i] = ONE_UP
        elif not tree[tree_dict[treeclass[i]]].is_leaf():
            print(tree.is_branch(tree_dict[treeclass[i]]),tree_dict[flatclass[i]]) # not sure why this needed
            if tree_dict[flatclass[i]] in tree.is_branch(tree_dict[treeclass[i]]): # flat class is one level down from tree class
                score[i] = ONE_DOWN
        # can imagine traversing more levels but prob not needed here
        # TODO: deal with SIMBAD 'Candidate' classes
        else:
            score[i] = NON_MATCH
        print(flatclass[i],treeclass[i],score[i])

    return(score)

def tree_read(infile):
    '''read SIMBAD class structure and put it into the treelib tree structure
       returns both a tree and the dictionary needed to translate btw nodeIDs and tags'''
    new_tree = Tree()
    tree_dict = {}
    dat = Table.read(infile, format='csv')
    numeric_id = dat['id']
    tag = dat['tag']
    for i, num_id in enumerate(numeric_id):
#        id_str = num_id.replace(".","")
        new_tree = add_to_tree(num_id, tag[i], new_tree)
        tree_dict[tag[i]] = num_id
    return(new_tree, tree_dict)

def get_parent(str_id):
    '''find the first occurrence of 00 in string, 
        replace previous 2 and all following digits with zeros'''
    double_zero = str_id.find("00")
    if double_zero == -1: # not found
        parent = str_id[:-3] 
    elif double_zero<4:
        parent = '00'
    else:
        parent = str_id[:double_zero-4]

    # pad out to length 11
    while len(parent) < 11:
        parent = parent+'.00'
    return(parent)

def add_to_tree(id, tag, tree):
    '''add an entry to treelib tree, using ab.cd.ef.gh as identifier
       where the parent of ab.cd.ef.gh is ab.cd.eg.00'''
    parent = get_parent(id)
    if parent == '00.00.00.00':
        tree.create_node(tag,id,parent=tree.root) # add to root
    elif tree.contains(parent): # add as child to parent
        tree.create_node(tag,id,parent=parent)
    else: # shouldn't get here, SIMBAD tree is in order
        print('create parent first')
    return(tree)
