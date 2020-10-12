"""takes matching files from nuc/ and pep/ directories, and if necessary
(e.g. nuc is transcript not CDS) finds the CDS to produce the provided pep sequence"""

from Bio import SeqIO, Seq, SeqRecord
import re
import os
import glob


def find_orfs(seq):
    seq_str = str(seq.seq)
    match_gen = re.finditer('ATG', seq_str)
    for match in match_gen:
        start = match.start()
        nostop = True
        out = ''
        i = start
        while nostop:
            codon = seq_str[i:(i + 3)]
            if codon in ['TAA', 'TGA', 'TAG' ]:
                nostop = False
            elif len(codon) < 3:
                break
            out += codon
            i += 3
        yield out


def clean_nucpep(nuc, pep, fixnuc):
    nucs = {}
    for record in SeqIO.parse(nuc, "fasta"):
        nucs[record.id] = record
    peps = {}
    for record in SeqIO.parse(pep, "fasta"):
        peps[record.id] = record
    
    for key in sorted(nucs.keys()):
        pep_provided = peps[key].seq
        orfs = list(find_orfs(nucs[key]))
        orfs = sorted(orfs, key=lambda x: -len(x)) # sort descending by length
        for orf in orfs:
            orfseq = Seq.Seq(orf)
            if pep_provided == Seq.Seq(orf).translate():
                found = orf
                break
            elif pep_provided == Seq.Seq(orf[:-3]).translate():
                found = orf[:-3]
                break
        nucs[key].seq = Seq.Seq(found)
    
    SeqIO.write(nucs.values(), fixnuc, "fasta")


def main():
    if not os.path.exists('fixnuc'):
        os.mkdir('fixnuc')
    
    nucs = os.listdir('nuc')
    for nuc in nucs:
        basef = re.sub('\.nuc', '', nuc)
        nuc_path = 'nuc/{}'.format(nuc)
        pep_path = glob.glob('pep/{}.*'.format(basef))[0]
        print(pep_path)
        fix_path = 'fixnuc/{}'.format(nuc)
    
        clean_nucpep(nuc_path, pep_path, fix_path) 

if __name__ == "__main__":
    main()
