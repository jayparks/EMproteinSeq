# EMproteinSeq
Running EM algorithm to group protein sequences

# Description
Grouping protein sequences by running EM algorithm.

# Usage instruction
* --input,     -i : input file containing protein sequences
* --seqlen,    -l : Length of sequences to use, (default: length of actual input sequence)
* --numfamily, -nf: Number of groups or families (default: 5)
* --iters,     -ii: Number of iterations to run (default: 100)
* --savefreq,  -sf: Save frequency of algorithm results (default: 10)

# Running demos
Sample running script: outputs class_label / softmembership
* python run_em.py -i subseqs10.txt
