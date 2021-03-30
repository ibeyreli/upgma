# UPGMA
Aim: Given n DNA sequences, 2 < n <25, in FASTA file, implement Guide
Tree Construction using the UPGMA algorithm. Output will be in Newick tree format
(http://evolution.genetics.washington.edu/phylip/newicktree.html).

# buildUPGMA.py

buildUPGMA.py is a program that performs UPGMA algorithm and outputs a guidance tree for given sequences.

## Installation

buildUPGMA.py uses Python 3 standard libraries.
If you have installed Python 3 before, no additional package is needed.
If not, you may install Python 3 from "https://www.python.org/downloads/".
Make sure you have selected the appropriate distribution for your operating system.

## Usage

You can run the program using terminal.

The syntax:

"python buildUPGMA.py --fasta sequences.fasta --gap ${gap_penalty} --match ${match_score}
 --mismatch ${mismatch_penalty} --out sequences.tree"
 
--fasta					: Indicator that the next input is the input file which
						contains the sequence to be considered during tree construction
sequences.fasta			: Sequence fasta file to be considered during tree construction
--match 				: Indicator that the next input is matching score
--gapop					: Indicator that the next input is gap opening penalty
--gapext				: Indicator that the next input is gap extention penalty, 
						usually less than gap opening penalty
--mismatch				: Indicator that the next input is mismatch penalty score
--out 					: Indicator that the next input is the name of the output 
						file which will contain the tree in Newick tree format
sequences.tree			: Output tree file in Newick tree format

## Examples:

sequences.fasta :
>A
CTAGATAATTGCCAGATGATCAAATTTATAT
>B
CTAGATAATCATGCTAGCTAGTGCACAAATTTATAT
>C
CTAGATAATTGGAATGTCGATCGATCG

Input >>> buildUPGMA --fasta sequences.fasta --match 5 --mismatch -3 --gapopen -8 --gapext -1 --out sequences.tree

Output >>>
sequences.tree :
(A:4.5,B:4.5):4.75,C:9.25);

sequences.fasta :
>A
CTAGA
>B
CTAGAT
>C
CTAGAA
>D
CTAGAAG

Input >>> buildUPGMA --fasta sequences.fasta --match 5 --mismatch -3 --gapopen -8 --gapext -1 --out sequences.tree

Output >>>
sequences.tree :
(A:0.5,B:0.5):0.25,(C:0.0,D:0.0):0.25);
