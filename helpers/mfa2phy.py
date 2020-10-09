__author__ = 'jmass'

import sys
from .fastahelper import FastaParser


def mfa2phy(mfa_in=None, phy_out=None):
    num_seq = 0
    alignment_length = None
    dct = {}
    for h, s in FastaParser().read_fasta(mfa_in):
        num_seq += 1
        if alignment_length is None:
            alignment_length = len(s.strip())
        else:
            if alignment_length != len(s.strip()):
                sys.stderr.write("Alignment lengths differ in file {} at {}!\n".format(mfa_in, h))
                sys.exit(1)
        dct[h] = s
    # write phylip-esque output
    with open(phy_out, 'w') as phy:
        phy.write("{}    {}\n".format(num_seq, alignment_length)) #dont ask... there were 4 spaces in the examples, that's why.
        for k, v in dct.items():
            offset = (30 - len(k))*" "
            phy.write("{}{}{}\n".format(k, offset, v))


def main():
    mfa2phy(sys.argv[1], sys.argv[1]+".phy")

if __name__=="__main__":
    main()
