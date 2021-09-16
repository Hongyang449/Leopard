import pyBigWig
import numpy as np
import os
import sys

def seq_to_hot(seq):
    import numpy as np
    import re
    seq=re.sub('B','N',seq)
    seq=re.sub('[D-F]','N',seq)
    seq=re.sub('[H-S]','N',seq)
    seq=re.sub('[U-Z]','N',seq)
    seq=seq.replace('a','A')
    seq=seq.replace('c','C')
    seq=seq.replace('g','G')
    seq=seq.replace('t','T')
    seq=seq.replace('n','N')
    Aseq=seq
    Aseq=Aseq.replace('A','1')
    Aseq=Aseq.replace('C','0')
    Aseq=Aseq.replace('G','0')
    Aseq=Aseq.replace('T','0')
    Aseq=Aseq.replace('N','0')
    Aseq=np.asarray(list(Aseq),dtype='float32')
    Cseq=seq
    Cseq=Cseq.replace('A','0')
    Cseq=Cseq.replace('C','1')
    Cseq=Cseq.replace('G','0')
    Cseq=Cseq.replace('T','0')
    Cseq=Cseq.replace('N','0')
    Cseq=np.asarray(list(Cseq),dtype='float32')
    Gseq=seq
    Gseq=Gseq.replace('A','0')
    Gseq=Gseq.replace('C','0')
    Gseq=Gseq.replace('G','1')
    Gseq=Gseq.replace('T','0')
    Gseq=Gseq.replace('N','0')
    Gseq=np.asarray(list(Gseq),dtype='float32')
    Tseq=seq
    Tseq=Tseq.replace('A','0')
    Tseq=Tseq.replace('C','0')
    Tseq=Tseq.replace('G','0')
    Tseq=Tseq.replace('T','1')
    Tseq=Tseq.replace('N','0')
    Tseq=np.asarray(list(Tseq),dtype='float32')
    hot=np.vstack((Aseq,Cseq,Gseq,Tseq))
    return hot 

# GRCh37
#chr_all=['chr' + str(i) for i in range(1,23)] + ['chrX']
chr_all=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
num_bp=np.array([249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560])
chr_len=dict(zip(chr_all,num_bp.tolist()))

path1='./grch37/'
path2='./dna_npy/'
os.system('mkdir -p ' + path2)
path3='./dna_bigwig/'
os.system('mkdir -p ' + path3)

for the_chr in chr_all:
    print(the_chr)
    f=open(path1 + the_chr + '.fa')
    f.readline()
    seq = ''
    for line in f:
        seq += line.rstrip()
    f.close()
    print(len(seq))
    seq_onehot=seq_to_hot(seq)
    np.save(path2 + the_chr, seq_onehot)

base_all=['A','C','G','T']
for i in np.arange(len(base_all)):
    print(base_all[i])
    bw = pyBigWig.open(path3 + base_all[i] + '.bigwig', 'w')
    bw.addHeader(list(zip(chr_all , num_bp)), maxZooms=0) # zip two turples
    for the_chr in chr_all:
        print(the_chr)
        x=np.load(path2 + the_chr + '.npy')
        # pad two zeroes
        z=np.concatenate(([0],x[i,:],[0]))
        # find boundary
        tmp1=np.where(np.diff(z)==1)[0]
        tmp2=np.where(np.diff(z)==-1)[0]
        starts=np.concatenate((tmp1, tmp2))
        starts.sort()
        ends=starts[1:]
        starts=starts[:-1]
        vals=np.zeros(len(starts))
        vals[np.arange(0,len(vals),2)]=1 # assume start with 0
        if starts[0]!=0: # if start with 1
            ends=np.concatenate(([starts[0]],ends))
            starts=np.concatenate(([0],starts))
            vals=np.concatenate(([0],vals))
        if ends[-1]!=chr_len[the_chr]: # if end with 0
            starts=np.concatenate((starts,[ends[-1]]))
            ends=np.concatenate((ends,[chr_len[the_chr]]))
            vals=np.concatenate((vals,[0]))
        # write
        chroms = np.array([the_chr] * len(vals))
        bw.addEntries(chroms, starts, ends=ends, values=vals)
    bw.close()


