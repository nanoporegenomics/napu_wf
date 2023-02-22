import numpy
import random
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


_nuc = numpy.array(["A", "T", "C", "G"])


def randSeq(length):
    seqArray = _nuc[[random.randint(0, 3) for i in range(length)]]
    return(MutableSeq("".join(seqArray)))


# parameters
N = 100000
snp_rate = .01
recomb_rate = .001
N_ont_reads = 100
ont_read_size = 11000

# "ancestral" sequence
anc = randSeq(N)
# write as "reference"
SeqIO.write([SeqRecord(anc, id='chr1')], "ref.fa", "fasta")

# 2 sets of parent haplotypes with some SNPs
parents = [[], []]
for bb in range(2):
    for ii in range(2):
        par = MutableSeq(str(anc))
        for pos in range(len(par)):
            if random.random() < snp_rate:
                new_base = _nuc[random.randint(0, 3)]
                while par[pos] == new_base:
                    new_base = _nuc[random.randint(0, 3)]
                par[pos] = new_base
        parents[bb].append(par)

# recombine parents into child haplotype pair
haps = ['', '']
for bb in range(2):
    cur_par = 0
    for pos in range(len(parents[bb][0])):
        haps[bb] += parents[bb][cur_par][pos]
        if random.random() < recomb_rate:
            cur_par = 1 - cur_par

# simulate ONT reads
ont_f = open('test_ont.fastq', 'wt')
readid = 0
for rr in range(N_ont_reads):
    hap = random.randint(0, 1)
    pos = random.randint(0, N - ont_read_size)
    read = haps[hap][pos:(pos+ont_read_size)]
    ont_f.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*ont_read_size))
    readid += 1
ont_f.close()
