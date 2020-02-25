# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: alindner1

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) reprfamiesented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    else:
        return 'C'

    pass

""" For this implementation, I created for loops to indicate if the value passed into the function (nucleotide)
is A, T, G, or C then return the complementary nucleotide (given above)

"""

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    complementstring = ''
    for i in dna:
        complementstring = complementstring + get_complement(i)
    final = complementstring[::-1]
    return final


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
            stop codon (TAG, TAA, TGA)
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATCCC")
    'ATGAGATCCC'
    """

    i = 0
    while i < len(dna):
        if dna[i:i+3] == 'TAG':
            return dna[:i]

        elif dna[i:i+3] == 'TAA':
            return dna[:i]


        elif dna[i:i+3] == 'TGA':
            return dna[:i]
        i = i + 3
    return dna

"""
reallyendbefore = 0
i = 0
while i < len(dna):
    if dna[i:i+3] == "TAG":
        endbefore = dna.index("TAG")
        reallyendbefore = endbefore +3
        return dna[:reallyendbefore]
    elif dna[i:i+3] == "TAA":
        endbefore = dna.index("TAA")
        reallyendbefore = endbefore +3
        return dna[:reallyendbefore]
    elif dna[i:i+3] == "TGA":
        endbefore = dna.index("TGA")
        reallyendbefore = endbefore + 3
        return dna[:reallyendbefore]
    else:
        return dna
        """
"""    for i in dna:
        if ((i == 'T') and ((i+1) == 'A') and ((i+3)== 'G')) or ((i == 'T') and ((i+1) == 'A') and ((i+3)== 'A'))  or ((i == 'T') and ((i+1) == 'G') and ((i+3)== 'A')):
            return toprint
        else:
            toprint = toprint + dna[i:i+3]
"""

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    # TODO: implement Sthis
    nonnested = []
    i = 0
    nonnested.append(rest_of_ORF(dna))
    while i < len(dna):
        #if dna[i:i+3] == 'ATG':
        #    i = 3
        orflength = len(rest_of_ORF(dna))

        i = orflength
    return nonnested
"""
    nonnested = []
    i = 0
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            orf = rest_of_ORF(dna[i::])
            nonnested.append(orf)
            i = i + len(orf)
        i = i + 3
    return nonnested

"""
        if "TAG" in dna:
            nonnested.append(rest_of_ORF(dna))
            startindexofstopcodon = dna.index("TAG")
            realstartindex = startindexofstopcodon + 3
            newstring = dna[realstartindex::]
            nonnested.append(rest_of_ORF(newstring))
        elif "TAA" in dna:
            nonnested.append(rest_of_ORF(dna))
            startindexofstopcodon = dna.index("TAA")
            realstartindex = startindexofstopcodon + 3
            newstring = dna[realstartindex::]
            nonnested.append(rest_of_ORF(newstring))
        elif "TGA" in dna:
            nonnested.append(rest_of_ORF(dna))
            startindexofstopcodon = dna.index("TGA")
            realstartindex = startindexofstopcodon + 3
            newstring = dna[realstartindex::]
            nonnested.append(rest_of_ORF(newstring))
        else:
            nonnested.append(rest_of_ORF(dna))

    return nonnested
"""


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    list2 = []

    list2.extend(find_all_ORFs_oneframe(dna[0:]))
    list2.extend(find_all_ORFs_oneframe(dna[1:]))
    list2.extend(find_all_ORFs_oneframe(dna[2:]))

    return list2


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    list3 = []
    list3.extend(find_all_ORFs(dna))
    list3.extend(find_all_ORFs(get_reverse_complement(dna)))
    return list3


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    list4 = find_all_ORFs_both_strands(dna)
    longer = ""
    for i in list4:
        if len(i) > len(longer):
            longer = i
    return longer



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        >>> longest_ORF_noncoding('ATGCCC', 10000)
        6
        """
        #shuffle_string find longest orf in it in, do it num trials number of times put in list and then longest
    i = 0
    longestones = []
    while num_trials > i:
        longestones.append((longest_ORF(shuffle_string(dna))))
        i = i + 1
    long = len(max(longestones))
    return long



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    protein = ''
    i = 0
    while i <= (len(dna) -3):
        protein += aa_table[dna[i:i+3]]
        i = i +3
    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        find all orfs on both strands
        remove all orfs less than threshold
        convert orfs into aa sequences
        return list
    """
    # TODO: implement this
    #orf on both strands remover less turn those into aa
    threshold = longest_ORF_noncoding(dna, 1500)
    genelist = []
    finallist = []
    genelist = find_all_ORFs_both_strands(dna)
    for i in genelist:
        if len(i) >= threshold:
            strand = coding_strand_to_AA(i)
            finallist.append(strand)
    return finallist

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(gene_finder, globals(), verbose = True)
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
