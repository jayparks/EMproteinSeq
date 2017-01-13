# EMproteinSeq
Running EM algorithm to group protein sequences

# Description
Grouping protein sequences by running EM algorithm.

# Usage instruction
Sample input running script : <code> python run_em.py -i subseqs10.txt </code>
* --input,     -i : input file containing protein sequences
* --seqlen,    -l : Length of sequences to use, (default: length of actual input sequence)
* --numfamily, -nf: Number of groups or families (default: 5)
* --iters,     -ii: Number of iterations to run (default: 100)
* --savefreq,  -sf: Save frequency of algorithm results (default: 10)

# Result
Sample output file : <file> subseqs10_iter0.txt </file>
* class_label / soft_membership
