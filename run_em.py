'''
Run EM algorithm to group protein sequences
'''

import sys
import argparse

import numpy as np


def normalize(array):
    return array / sum(array)

def create_data(input_path, seqlen):
    fin = open(input_path, 'r')
    lines = fin.readlines()
    # exclude first line which is for sequence detail
    lines = lines[1:]
    if seqlen <= 0:
        # if seqlen <= 0, use actual sequence length (excluding '\n')
        seqlen = len(lines[0])-1
    array = np.array([list(line.strip()[:seqlen]) for line in lines if line.strip() != ''])
    print("created data with size (%d, %d)" % (array.shape[0], array.shape[1]))
    
    return array
        
def create_dictionary(data):
    dict = {}
    idx = 0
    for char in np.nditer(data):
        if char.tolist() not in dict:
            dict[char.tolist()] = idx
            idx += 1
    print("created dictionary with %d characters" % len(dict))

    return dict

def create_membership(seqsize, numfamily):
    
    print("created membership matrix of size (%d, %d)" % (seqsize, numfamily))
    return np.apply_along_axis(
        normalize, axis=1, arr=np.random.rand(seqsize, numfamily))

def create_prob(numfamily, numchar, seqlen):
    
    print("created probability matrix of size (%d, %d, %d)" % (numfamily, numchar, seqlen))
    return np.apply_along_axis(
        normalize, axis=1, arr=np.random.rand(numfamily, numchar, seqlen))

def compute_membership(seq, family, dict, prob_mat):
    likelihood = 1.0
    for col, char in enumerate(seq):
        likelihood *= prob_mat[family][dict[char]][col]
    return likelihood
    
def save_result(input_path, idx, membership):
    output_path = input_path[:input_path.index(".")]
    f = open("%s_iter%d.txt" % (output_path, idx), 'w')
    for seq in membership:
        f.write("%d %.2f\n" % (np.argmax(seq), seq[np.argmax(seq)]))
    

def main(input_path, seqlen=0, numfamily=5, iters=100, savefreq=10):
    
    # build data array for whole sequence
    data = create_data(input_path, seqlen)

    seqlen = data.shape[1]   
 
    # build dictionary for amino acid from input data
    dict = create_dictionary(data)
    
    # initialize membership array for entire sequences: seqsize x numfamily
    membership = create_membership(len(data), numfamily)
    
    # initialize probability matrix for each family: numfamily x len(dict) x seqlen
    prob_mat = create_prob(numfamily, len(dict), seqlen)
    
    print("Running EM algorithm..")
    for idx in xrange(iters):
        if idx % 10 == 0:
            print("  %dth iter.." % idx)
        if idx % savefreq == 0:
            print("  saving results on %dth iter.." % idx)
            save_result(input_path, idx, membership)
        
        # Expectation
        for family in xrange(numfamily):
            for col in xrange(seqlen):
                for char in dict:
                    mask = data[:, col] != char
                    # if the $col of data contains no character $char,
                    if np.ma.all(mask):
                        prob_mat[family][dict[char]][col] = 0
                    # if the column of data contains at least one $char
                    else:
                        masked_sum = np.ma.masked_array(membership[:, family], mask).sum()
                        prob_mat[family][dict[char]][col] = masked_sum / sum(membership[:, family])
        # Maximization
        for seqidx, seq in enumerate(data):
            for family in xrange(numfamily):
                membership[seqidx, family] = compute_membership(seq, family, dict, prob_mat)
            membership[seqidx] = np.apply_along_axis(normalize, axis=0, arr=membership[seqidx])
            
        
    print("DONE..")
    


if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, help="Input file path")
    parser.add_argument('--seqlen', '-l', type=int, default=0,
                        help="Length of sequences to use, (default: length of actual input sequence")
    parser.add_argument('--numfamily', '-nf', type=int, default=5,
                        help="Number of families")
    parser.add_argument('--iters', '-ii', type=int, default=100, 
                        help="Number of iterations to run")
    parser.add_argument('--savefreq', '-sf', type=int, default=10,
                        help="Save frequency of algorithm results")
    
    args = parser.parse_args()

    main(args.input, args.seqlen, args.numfamily, args.iters, args.savefreq)
    
    
