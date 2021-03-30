# -*- coding: utf-8 -*-
"""
CS481 Fall 2019 Homework 5
UPGMA
Source Code File
Created on Wed Dec 11 22:42:40 2019
@author: İlayda Beyreli 201801130

Notification:
    This Python file contains function descriptions and/or algorithms
    that have been coded by the author for previous asignments 
    with or withour some modifications.
"""
import sys # To get arguments from terminal, included in the standard library

def readfa_multi(file_name):
    # Reading .FA files into lists
    # Special setup for given assignment
    seqs  = []
    names = []
    with open(file_name, "r") as fin:
        data =fin.read().splitlines()
        while data:
            seqs.append(data[1])
            names.append(data[0][1:])
            data = data[2:]
    return names, seqs
#=============================================================================

def Newick(clusters, ntree=""):
    # Recurrent Newick Tree Formatter
    while isinstance(clusters,list):
        h = clusters[2]/2
        if isinstance(clusters[0], list):
            ntree += Newick(clusters[0])+":"+str(h-clusters[0][2]/2)+","
        else:
            ntree += "("+clusters[0]+":"+str(h)+","
        if isinstance(clusters[1], list):
            ntree += Newick(clusters[1])+":"+str(h-clusters[0][2]/2)+")"
        else:
            ntree += clusters[1]+":"+str(h)+")"
        clusters = clusters[2]
    return ntree

# =============================================================================

def zeros(r,c):# Matrix of zeros as list of lists 
    return [[0 for i in range(c)] for i in range(r)]

def empty(r,c):# Matrix of negative inf as list of lists 
    return [[float("inf") for i in range(c)] for i in range(r)]

def matprint(matrix): # Print function for matrices stored as list of lists
    for i in range(len(matrix)):
            print(matrix[i],"\n")
    return 0

def argmin(l):
    return l.index(min(l))

def argmax(l):
    return l.index(max(l))

# =============================================================================
    
def naive_score(c1,c2,match=0,gapop=1,mismatch=2):
    # Scoring Function
    # defaults are given according to edit distance formulation
    if c1 == c2:
        return match
    elif c1 =='-' or c2 =='-':
        return gapop
    else:
        return mismatch

def edit_distance(seq1,seq2):
    # Initialize nm table
    n = len(seq1)
    m = len(seq2)
    S = zeros(n+1,m+1)
    for i in range(n+1):
        S[i][0] = i
    for i in range(m+1):
        S[0][i] = i
    # Fill the scoring table 
    for i in range(1,n+1):
        for j in range(1,m+1):
            match = S[i-1][j-1]+naive_score(seq1[i-1],seq2[j-1])
            del1 = S[i-1][j]+naive_score(seq1[i-1],'-')
            in1 = S[i][j-1]+naive_score('-',seq2[j-1])
            S[i][j] = min(match,del1,in1)
    return S[-1][-1]

def direction(a1,a2,a3,a4):
    switcher={0:'d',
            1:'u',#del
            2:'l',#in
            3:'d'}
    l = [a1,a2,a3,a4]
    m = max(a1,a2,a3,a4)
    m = l.index(m)
    return switcher.get(m)

def affine_dist(seq1,seq2,match,gapop,gapext,mismatch):
    # Needleman-Wunsh Affine Alignment Algorithm and Edit Distance
    # Create tables
    n = len(seq1)
    m = len(seq2)  
    G = zeros(n+1,m+1)
    E = zeros(n+1,m+1)
    F = zeros(n+1,m+1)
    V = zeros(n+1,m+1)
    L = zeros(n+1,m+1)
    # Initialize according to global alignment    
    l = -float("inf")
    G[0][0] = -float("inf")
    for j in range(m+1):
        F[0][j] = -float("inf")
    for i in range(n+1):
        E[i][0] = -float("inf")
    for j in range(1,m+1):
        V[0][j] = gapop+j*gapext
    for i in range(1,n+1):
        V[i][0] = gapop+i*gapext
    # Fill the DP table   
    for i in range(1,n+1):
        for j in range(1,m+1):
            a = E[i][j-1]+gapext
            b = G[i][j-1]+gapop+gapext
            c = F[i][j-1]+gapop+gapext
            E[i][j]=max(a,b,c)
            
            a = F[i-1][j]+gapext
            b = G[i-1][j]+gapop+gapext
            c = E[i-1][j]+gapop+gapext
            F[i][j]=max(a,b,c)
            
            G[i][j] = V[i-1][j-1]+naive_score(seq1[i-1],seq2[j-1],match,gapop,mismatch)
            V[i][j] = max(G[i][j],E[i][j],F[i][j])
            L[i][j] = direction(G[i][j],E[i][j],F[i][j],l)

    # Traversing
    aseq1=[]
    aseq2=[]
    i = n
    j = m
    while j > -1 and i > -1:
        # print(L[i][j]) # Debugging
        if L[i][j] == 'd':
            aseq1.append(seq1[i-1])
            aseq2.append(seq2[j-1])
            i-=1
            j-=1
        elif L[i][j] == 'l':
            aseq1.append('-')
            aseq2.append(seq2[j-1])
            i-=1
        elif L[i][j] == 'u':
            aseq1.append(seq1[i-1])
            aseq2.append('-')    
            j-=1
        else:
            i-=1
            j-=1
    # print(aseq1,aseq2) # Debugging
    aseq1.reverse()
    aseq1=''.join(aseq1)
    aseq2.reverse()    
    aseq2=''.join(aseq2) 
    # Return the edit distance between pairwise-alignments
    return edit_distance(aseq1,aseq2)

# =============================================================================

def upgma(seqs, names, match, gapop, gapext, mismatch):
    clusters  = names
    n = len(seqs)
    # Creating the distance matrix
    dists = empty(n,n)
    # Computing the distance matrix.
    for i in range(n):
        for j in range(i):
            dists[i][j]=affine_dist(seqs[i],seqs[j],match,gapop,gapext,mismatch)
    # While there are at least two clusters, do following...
    while len(dists) > 1:
        # get the closest distance
        minscore = min([min(row) for row in dists])
        score_r = argmin([min(row) for row in dists])
        score_c = argmin(dists[score_r])
        # print(minscore, score_r,score_c)
        # get nodes connected by shortest edge
        c1 = clusters[score_c]
        c2 = clusters[score_r]
        clusters.remove(c1)
        clusters.remove(c2)
        # insert merged cluster
        clusters.insert(min(score_r,score_c),[c1,c2,minscore])
        # recalculate all other distances
        dists.pop(max(score_r,score_c))
        for row in dists:
            row[min(score_r,score_c)]=(row[min(score_r,score_c)]+row[max(score_r,score_c)])/(len(c1)+len(c2))
            row.pop(max(score_r,score_c))
        #print(c1,c2)
        #print(clusters)
        #matprint(dists)
    return clusters[0]

# =============================================================================
# Start the main program
# =============================================================================
if __name__ == '__main__':
    
    fasta_file = [sys.argv[i] for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--fasta"][0]
    out = [sys.argv[i] for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--out"][0]
    gapop = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--gapopen"][0]
    gapext = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--gapext"][0]
    match = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--match"][0]
    mismatch = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--mismatch"][0]
    
    """
    # Debugging :p
    fasta_file = "hw5dummy.fa"
    out = "hw5example.txt"
    gapop = -8
    gapext = -1
    match = 5
    mismatch = -3
    """
    names, seqs = readfa_multi(fasta_file)    
    t = upgma(seqs, names, match,gapop,gapext, mismatch)
    # print(t) # Debugging
    ntree = Newick(t)+";"
    # print(ntree) # Debugging
    with open(out,"+w") as fout:
        fout.write(ntree)