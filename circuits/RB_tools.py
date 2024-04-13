# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:02:02 2021

@author: Karl Mayer
"""

import numpy as np
import pickle

pi = np.pi

# load TQ Clifford group
TQ_Clifford_group = pickle.load(open('TQ_Clifford_group.p', 'rb'))
Clifford_group_list = list(TQ_Clifford_group.keys())


def generate_group(generators):
    """ generators (dict): str/np.array key/value pairs
                           ex. {('a'):np.array, ('b'):np.array, ...}
                           
        returns dictionary of full group
        
        Note: group alphabet must only contain length-1 letters
    """
    
    # initialize group
    group = dict(generators)
    
    # add identity
    d = len(list(generators.values())[0])
    group['I'] = np.identity(d)
    
    # look for new group elements
    finished = False
    L = 1 # length of group elements searched over so far
    
    while finished == False:
        
        new_g = {}
        for g_str in group:
            #if g_str != 'I' and len(g_str) == L:
            if g_str != 'I' and g_str.count('*') == L-1:
                g = group[g_str]
                for h_str in generators:
                    h = group[h_str]
                    #gh_str = g_str + h_str
                    gh_str = g_str + '*' + h_str
                    gh = g @ h
                    if is_in_group(gh, group) == False and is_in_group(gh, new_g) == False:
                        new_g[gh_str] = gh
        
        # add the new group elements
        if new_g == {}:
            finished = True
        else:
            for gh_str in new_g:
                group[gh_str] = new_g[gh_str]
                print(len(group), gh_str)
            L += 1
                    
    return group


def is_in_group(g, group):
    """ g: np.array
        group (dict): str/np.array key/value pairs
    """
    
    d = len(g)
    
    g_dg = np.conj(np.transpose(g))
    
    answer = False
    
    # search through group for elements equivalent to g
    for h_str in group:
        h = group[h_str]
        dist = 1 - (np.abs(np.trace(g_dg @ h))/d)**2
        if dist < 10**(-8):
            answer = True
    
    return answer


def make_RB_seq(seq_len):
    
    # sample random Clifford group elements
    group_elements = list(np.random.choice(Clifford_group_list, size=seq_len))

    # compute inverse gate
    circ_U = np.identity(4)
    for group_el in group_elements:
        g = TQ_Clifford_group[group_el]
        circ_U = g @ circ_U
        
    # search through group for elements equivalent to U_inv
    def find_inverse(U):
        #U_inv = np.conj(U.T)
        for group_el in Clifford_group_list:
            g = TQ_Clifford_group[group_el]
            dist = 1 - (np.abs(np.trace(U @ g))/4)**2
            if dist < 10**(-8):
                break

        return group_el

    g_inv = find_inverse(circ_U)
    
    sequence = group_elements + [g_inv]
    
    return sequence


