import numpy as np
import pysam
from pysam import VariantFile
import random
import sys, os

reference = pysam.FastaFile(str(argv.sys[1]))
samfile=pysam.AlignmentFile(str(argv.sys[2]), "rb")
samfile2=pysam.AlignmentFile(str(argv.sys[3]), "rb")
out=pysam.AlignmentFile("out.bam", "wb",template= samfile2)
bcf_in = VariantFile(str(argv.sys[4]))


def prob_matrix(samfile, reference, num_of_pos):
    positions = np.arange(num_of_pos+1)
    probability, nucl = np.zeros((num_of_pos+1,4,4)), np.zeros((num_of_pos+1,4,4))

    for position in positions:
        for region in reference.references:
            nucleotides = []
            for nucleotide in reference.fetch(region):
                nucleotides.append(nucleotide)

            indexesA = [i for i, x in enumerate(nucleotides) if x == 'A']
            indexesT = [i for i, x in enumerate(nucleotides) if x == 'T']
            indexesG = [i for i, x in enumerate(nucleotides) if x == 'C']
            indexesC = [i for i, x in enumerate(nucleotides) if x == 'G']
            all_index = [indexesA,indexesT,indexesG,indexesC] 

            for j in range(4):# for each ref nucl ATGC
                for pileupcolumn in samfile.pileup(region):
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip and pileupcolumn.pos in all_index[j]: 
                            if pileupread.query_position == position:
                                if pileupread.alignment.query_sequence[position] == "A":
                                    nucl[position][j][0]+=1
                                elif pileupread.alignment.query_sequence[position] == "T":
                                    nucl[position][j][1]+=1
                                elif pileupread.alignment.query_sequence[position] == "C":
                                    nucl[position][j][2]+=1
                                elif pileupread.alignment.query_sequence[position] == "G":
                                    nucl[position][j][3]+=1
    for position in range(num_of_pos+1):
        for i in range(4):
            for j in range(4):
                probability[position][i][j]+=(nucl[position][i][j]/sum(nucl[position][i]))
    return probability


matrix = prob_matrix(samfile,reference, 39)


def get_nucl(matrix,i,j):
    if matrix[i][j].argmax() in (3,7,11,15):
        return 'G'
    elif matrix[i][j].argmax() in (2,6,10,14):
        return 'C'
    elif matrix[i][j].argmax() in (1,5,9,13):
        return 'T'
    else:
        return 'A'


def clear_error(samfile2, reference):
    read_list=[]
    for region in reference.references:
        nucleotides = []
        for nucleotide in reference.fetch(region):
            nucleotides.append(nucleotide)

        indexesA = [i for i, x in enumerate(nucleotides) if x == 'A']
        indexesT = [i for i, x in enumerate(nucleotides) if x == 'T']
        indexesG = [i for i, x in enumerate(nucleotides) if x == 'C']
        indexesC = [i for i, x in enumerate(nucleotides) if x == 'G']

        variants=[]
        positions=[]
        for rec in bcf_in.fetch(region):
            variants.append(rec.alleles)
            positions.append(rec.pos)

        for read in samfile2.fetch(region):
            strread=read.query_sequence
            for i in range(len(strread)):
                ref_pos=read.pos+i
                if ref_pos in positions and strread[i] in variants[positions.index(ref_pos)]:
                    strread_list=list(strread)
                    strread_list[i] = random.choice(variants[positions.index(ref_pos)])
                    strread="".join(strread_list)
                if ref_pos in indexesA:
                    strread_list=list(strread)
                    strread_list[i]=get_nucl(matrix,i,0)
                    strread="".join(strread_list)
                elif ref_pos in indexesC:
                    strread_list=list(strread)
                    strread_list[i]=get_nucl(matrix,i,1)
                    strread="".join(strread_list)
                elif ref_pos in indexesT:
                    strread_list=list(strread)
                    strread_list[i]=get_nucl(matrix,i,2)
                    strread="".join(strread_list)
                else:
                    strread_list=list(strread)
                    strread_list[i]=get_nucl(matrix,i,3)
                    strread="".join(strread_list)
            read_list.append(strread)
    return read_list


reads = clear_error(samfile2, reference)

samfile2=pysam.AlignmentFile(str(argv.sys[3])), "rb")
out=pysam.AlignmentFile("out.bam", "wb",template= samfile2)
i=0
for s in samfile2:
    s.query_sequence = reads[i]
    i+=1
    out.write(s)
out.close()
