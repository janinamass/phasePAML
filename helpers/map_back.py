#########################
# last update:
# Fri 11 Apr 2014 12:21:41 PM CEST
# [JMass]
#########################
import getopt
import imp
import sys
"""
Maps coordinates eg from an ungapped alignment back to their original positions.
Needs:
gapped alignment (for now: paml format),
coordinate from ungapped alignment to map back relative to alignment.
"""
def usage():
    print ("""
    #############################################
    # example usage:
    # python phasePaml_mapBack.py -a alignment.gapped.paml -p 1,2,3
    ############################################

    general options:
    -a, --gapped_alignment= aln.paml   paml format (eg from pal2nal)
    -p, --positions NUM,NUM,NUM         eg 15,16,17,21

    -h, --help                  prints this
    """)
    sys.exit(2)

class Alignment(object):
    def __init__(self, alignment, alignment_format="PAML", gapped=True):
        self._alignment = alignment
        self._alignment_format = alignment_format
        self._sequences = self.readAlignedSequences()
        self.gapped = gapped
        self.gapPositions = self.calcGapped()

    def readAlignedSequences(self):
        """read PAML format"""
        headerline = True
        header = ""
        isID = False
        isSeq = False
        curSeq = ""
        curHead = ""
        resdct ={}
        if self._alignment_format.upper() == "PAML":
            with open(self._alignment) as al:
                for l in al:
                    if headerline:
                        header = l.strip()
                        headerline = False
                        isID = True
                    elif isID:
                        isID=False
                        curHead = l.strip()
                        currSeq = ""
                        isSeq = True
                    elif isSeq:
                        if l.strip() == curHead or l.strip()=="":
                            pass
                        elif (len(curSeq)<int(header.split(" ")[-1])):
                            curSeq+=l.strip()
                        else:
                            isSeq=True
                            resdct[curHead] = curSeq
                            curSeq = ""
                            curHead = l.strip()
                resdct[curHead]=curSeq
                print(resdct)
        else:
            sys.stderr.write("PAML format expected.\n")
            system.exit(2)
        return(resdct.copy())

    def calcGapped(self):
        res = []
        if self.gapped:
            ks = self._sequences.keys()
            for pos in range(0,len(self._sequences[ks[0]])):
                if [self._sequences[x][pos] for x in ks if self._sequences[x][pos] =="-"] :
                    res.append(pos)
            return(res)
        else:
            self.gapPositions = None

    def getCoordinatesMap(self):
        gapped = {}
        ungapped = {}
        ks = self._sequences.keys()
        if self.gapped:
            for i in range(0,len(self._sequences[ks[0]])):
                    gapped[i] = self.getMappedPos(i)
                    ungapped[self.getMappedPos(i)] = i
            return(gapped.copy(), ungapped.copy())
        else:
            return(None)

    def getMappedPos(self, position):
        res = None
        if (position in self.gapPositions):
            return None
        else:
            return(position - len([x for x in self.gapPositions if x < position]))

    #sequence specific
    def getMappedPosRef(self,k):
        seq = self._sequences[k] #seq with gaps
        gaps = []
        mappedRef = {}
        #map gapped positions
        for e,i in enumerate(seq):
            if i=="-":
                gaps.append(e)
            else:
                mappedRef[e] = e-len(gaps)
        return(mappedRef.copy())


def main():
    gappedAlignment = None
    pos =[]
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "a:p:h", ["gapped_alignment=","positions=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-a", "--gapped_alignment"):
            gappedAlignment = a
        elif o in ("-p", "--positions"):
            pos = [int(i) for i in a.strip().split(",")]

        elif o in ("-h","--help"):
            usage()
        else:
            assert False, "unhandled option"
    ########################################
    cgap = None
    cungap = None
    if gappedAlignment:
        a = Alignment(alignment = gappedAlignment,gapped=True)
    for p in pos:
        try:
            cgap = a.getCoordinatesMap()[0][p]
            cungap = a.getCoordinatesMap()[1][p]
        except KeyError as k:
            print(k)

        try:
            #print("#gapped to ungapped\n{}\t{}".format(p,cgap))
            print("#id\tungapped\tgapped\toriginal")
            for k in a._sequences:
                orig = a.getMappedPosRef(k)[cungap]
                print("{}\t{}\t{}\t{}".format(k,p,cungap,orig))
                #print(a._sequences[k],a.getMappedPosRef(k)[cungap])
        except KeyError as e:
            sys.stderr.write(str(e))

if __name__ == "__main__":
    main()
