#!/usr/bin/env python
import sys
import re

from ete2 import Tree
from ete2 import EvolTree
from ete2 import TreeStyle
from ete2 import NodeStyle
from ete2 import TextFace
import os

def ctlTemplate(seqfile=None, treefile=None, outfile=None, runmode="0", \
                seqtype="1", model=None, fix_kappa="0", kappa="2", \
                fix_omega=None, omega=None, NSsites=None, CodonFreq="2"):
    seqfile = "seqfile = {}\n".format(seqfile)
    treefile = "treefile = {}\n".format(treefile)
    outfile = "outfile = {}\n".format(outfile)
    runmode = "runmode = {}\n".format(runmode)
    seqtype = "seqtype = {}\n".format(str(seqtype))
    CodonFreq = "CodonFreq = {}\n".format(CodonFreq)
    model = "model = {}\n".format(model)
    NSsites = "NSsites = {}\n".format(NSsites)
    fix_kappa = "fix_kappa = {}\n".format(str(fix_kappa))
    # 0 * 1: kappa fixed, 0: kappa to be estimated
    kappa = "kappa = {}\n".format(str(kappa))
    #2 * initial or fixed kappa
    fix_omega = "fix_omega = {}\n".format(fix_omega)  #* 1: fixed, 0: estimate
    omega = "omega = {}\n".format(str(omega))  #1 * initial or fixed omega, for codons or codon-based AAs

    runmodeExplanation = """
    *********************
    * 0: user tree;
    * 1: semi-automatic;
    * 2: automatic
    * 3: StepwiseAddition;
    * (4,5):PerturbationNNI;
    * -2: pairwise
    *********************
    """
    seqtypeExplanation = """
    * 1:codons; 2:AAs; 3:codons-->AAs
    """
    CodonFreqExplanation = """
    *0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    """
    modelExplanation = """
    **********************
    * models for codons:
    * 0:one,
    * 1:b (free ratio),
    * 2:2 or more dN/dS ratios for branches
    * When model = 2, you have to group branches on the tree into branch groups
    * using the symbols # or $ in the
    **********************
    """
    NSsitesExplanation = """
    **********************
    * 0:one w;1:neutral;
    * 2:selection; 3:discrete;
    * 4:freqs; 5:gamma;
    * 6:2gamma;7:beta;
    * 8:beta&w;9:beta&gamma;
    * 10:beta&gamma+1; 11:beta&normal>1;
    * 12:0&2normal>1; 13:3normal>0
    **********************
    """
    stuff = """
    *icode = 0 * 0:universal code; 1:mammalian mt; 2-11: ...
    *Mgene = 0 * 0:rates, 1:separate
    ***************************************
    *fix_alpha = 1 *0: estimate gamma shape parameter; 1: fix it at alpha
    *alpha = 0.1 *initial or fixed alpha, 0:infinity (constant rate)
    *Malpha = 0 * different alphas for genes
    *ncatG = 8 *# of categories in dG of NSsites models
    *clock = 0 * 0:no clock, 1:clock; 2:local clock; 3:TipDate
getSE = 0 * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0 * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .5e-6
cleandata = 1 * remove sites with ambiguity data (1:yes, 0:no)
method = 0 * 0: simultaneous; 1: one branch at a time
    **************
    **************

    """
    template = seqfile + treefile + outfile
    template += runmode + seqtype + CodonFreq + model
    template += NSsites + fix_kappa + kappa
    template += fix_omega + omega
    template += stuff

    return (template)


def generateCtl(model, treefile, seqfile, outfile, generateOther=False):
    treefile = os.path.basename(treefile)
    seqfile = os.path.basename(seqfile)
    resfile = os.path.basename(outfile)
    #outfile = os.path.basename(outfile)
    print(outfile, seqfile, treefile)
    print("GENERATE CTL", treefile, seqfile)
    if not outfile:
        outfile = treefile
    # ######## Basic model ########################
    if (model.upper() == "M0"):
        outsuff0 = ".M0_e05"
        outsuff1 = ".M0_e10"
        outsuff2 = ".M0_f10"

        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model="0", NSsites="0",
                           fix_omega="0", omega="0.5")
        ctl1 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff1, model="0", NSsites="0",
                           fix_omega="0", omega="1")
        ctl2 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff2, model="0", NSsites="0",
                           fix_omega="1", omega="1")
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)
        with open(outfile + outsuff1 + ".ctl", 'w') as out:
            out.write(ctl1)
        with open(outfile + outsuff2 + ".ctl", 'w') as out:
            out.write(ctl2)
    ########## /Basic model ######################

    ##########  Branch models#####################
    elif (model.upper() == "FR"):
        outsuff0 = ".FR_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model="1", NSsites="0",
                           fix_omega=0, omega=1)
        #free ratio model  - use is discouraged
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)

    elif (model.upper() == "BM"):
        outsuff0 = ".BM_f1"
        outsuff1 = ".BM_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model="2", NSsites="0",
                           fix_omega=1, omega=1)
        ctl1 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff1, model="2", NSsites="0",
                           fix_omega=0, omega=1)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)
        with open(outfile + outsuff1 + ".ctl", 'w') as out:
            out.write(ctl1)

    ############ /Branch models ######################

    ############  Site models   ######################
    # M1a
    elif (model.upper() == "M1A"):
        #M1aNearlyNeutral
        outsuff0 = ".M1a_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=0, NSsites=1,
                           fix_omega=0, omega=1)
        if generateOther:
            generateCtl("M2A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)

    #M2a
    elif (model.upper() == "M2A"):
        #M2a (PositiveSelection)
        outsuff0 = ".M2a_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=0, NSsites=2,
                           fix_omega=0, omega=1)
        if generateOther:
            generateCtl("M1A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)

    #M7
    elif (model.upper() == "M7"):
        #M7 (beta)
        outsuff0 = ".M7_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=0, NSsites=7,
                           fix_omega=0, omega=1)
        if generateOther:
            generateCtl("M8", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)
    #M8
    elif (model.upper() == "M8"):
        outsuff0 = ".M8_e1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=0, NSsites=8,
                           fix_omega=0, omega=1)
        if generateOther:
            generateCtl("M7", treefile, seqfile, outfile, False)
            generateCtl("M8A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)
    #M8a
    elif (model.upper() == "M8A"):
        outsuff0 = ".M8a"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=0, NSsites=8,
                           fix_omega=1, omega=1)
        if generateOther:
            generateCtl("M8", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)

    #BranchSite models
    #modified model A Null Hypothesis
    elif (model.upper() == "AH0"):
        outsuff0 = ".Ah0"
        #null
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=2, NSsites=2,
                           fix_omega=1, omega=1)
        if generateOther:
            generateCtl("AH1", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)
    #modified model A Alternative
    elif (model.upper() == "AH1"):
        outsuff0 = ".Ah1"
        ctl0 = ctlTemplate(treefile=treefile, seqfile=seqfile, outfile=resfile + outsuff0, model=2, NSsites=2,
                           fix_omega=0, omega=1.5)
        if generateOther:
            generateCtl("AH0", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile + outsuff0 + ".ctl", 'w') as out:
            out.write(ctl0)

            #PAML manual:
            #To calculate the p value based on this mixture distribution,
            #you calculate the p value using X_1^2, and then divide it by 2.
            #Suppose your 2DELTAl= 2.71, you will get 0.10 from X_1^2 ,
            #the the p value for the mixture is 0.10/2 = 5%.
            #We recommend that you use X_1^2 (with critical values 3.84 and 5.99)
            #instead of the mixture to guide against violations of model assumptions.


def fakeUnroot(tree, bgcolor="white"):
    """ As a convention, in newick format,
    an unrooted bifurcating tree
    is represented by a 'root' that is trifurcating,
    and all its other nodes bifurcating.
    This function only hides the root and the branch leading
    to it while drawing
    """
    nds = NodeStyle()
    nds["hz_line_color"] = "white"
    nds["size"] = 0
    rt = tree.get_tree_root()
    rt.set_style(nds)
    return tree


def read_tree(tree):
    t = Tree(tree)
    for node in t.traverse():
        if node.name.startswith("""'"""):
            node.name = node.name.replace("""'""", "")
            node.name = node.name.replace(" ", "_")
    return t.write(format=9)


def labelForPaml(unlabelledTreeAsString, listOfNodes, tree):
    print("ULT", unlabelledTreeAsString)
    t = EvolTree(unlabelledTreeAsString)
    marks = []
    count = 1
    for i in listOfNodes:
        marks.append("#" + str(count))
        count += 1
    t.mark_tree(listOfNodes, marks=marks)
    print(t.write(format=9))
    outfile = tree + "." + "_".join(listOfNodes)
    with open(outfile, 'w') as out:
        out.write(t.write())
    return (outfile)


def rgb(c):
    split = (c[0:2], c[2:4], c[4:6])
    return [int(x, 16) for x in split]


def label_regex(unlabeled_tree, regex, outfile, depth=4):
    pattern = re.compile(regex)
    t = EvolTree(unlabeled_tree)
    marks = []
    count = 1
    outfiles = []
    ts = TreeStyle()
    ts.mode = "c"
    nsMatch = NodeStyle()  # match
    nsMatch["fgcolor"] = "blue"
    nsMatch["size"] = 10
    nsBG = NodeStyle()
    nsBG["fgcolor"] = "black"
    nsBG["size"] = 0
    nsFG = []

    tolabelreg = []
    for i in range(0, depth):
        nsFG.append(NodeStyle())
        nsFG[i]["size"] = 10
        nsFG[i]["fgcolor"] = "blue"

    isroot = True
    for node in t.traverse():
        node.set_style(nsBG)
        if node.is_root():
            print("root")
            node.unroot()
            node._support = None

    for node in t.get_descendants():
        node.add_face(TextFace(node.node_id), column=0)

    # traverse and match
    for node in t.traverse():
        if re.match(pattern, node.name):
            node.set_style(nsMatch)
            n = node
            try:
                for i in range(0, depth):
                    n = n.up
                    n.set_style(nsFG[i])
                    marks.append("#" + str(count))
                    t.mark_tree([str(count)], marks=marks)
                    #just label everything with #1
                    print(t.write())

                    tolabelreg.append(str(n.node_id))

                    outfile = outfile + "." + "_".join([str(n.node_id)])
                    with open(outfile, 'w') as out:
                        out.write(t.write())
                    outfiles.append(outfile)

            except AttributeError:
                pass

            marks.append("#" + str(count))
            print(count)
            t.mark_tree([str(count)], marks=marks)
            #just label everything with #1
            print(t.write())

            outfile = tree + "." + "_".join([str(node.node_id)])
            with open(outfile, 'w') as out:
                out.write(t.write())
                outfiles.append(outfile)
    #t.show()
    t.render(tree + ".png", tree_style=ts)
    return (outfiles, tolabelreg)


def label_regex(unlabeled_tree, regex, treefile, depth=4, model_list=None, paml_msa=None, outfile=None):
    pattern = re.compile(regex)
    t = EvolTree(unlabeled_tree)
    marks = []
    count = 1
    outfiles = []
    ts = TreeStyle()
    ts.mode = "c"
    nsMatch = NodeStyle()  # match
    nsMatch["fgcolor"] = "blue"
    nsMatch["size"] = 8
    nsBG = NodeStyle()
    nsBG["fgcolor"] = "black"
    nsBG["size"] = 0
    nsFG = []

    tolabelreg = []
    for i in range(0, depth):
        nsFG.append(NodeStyle())
        nsFG[i]["size"] = 8
        nsFG[i]["fgcolor"] = "blue"

    #isroot = True
    #for node in t.traverse():
    #    node.set_style(nsBG)
    #    if node.is_root():
    #        print("root")
    #       node.unroot()
    #       node._support = None

    for node in t.get_descendants():
        node.add_face(TextFace(node.node_id), column=0)

    # traverse and match
    for node in t.traverse():
        if re.match(pattern, node.name):
            print("MATCH", node.name, node.node_id)
            node.set_style(nsMatch)
            n = node
            try:
                for i in range(0, depth):
                    n = n.up
                    n.set_style(nsFG[i])
                    marks.append("#" + str(count))
                    #print(count)
                    t.mark_tree([str(count)], marks=marks)
                    #just label everything with #1
                    tolabelreg.append(str(n.node_id))

                    outfile = treefile + "." + "_".join([str(n.node_id)])
                    with open(outfile, 'w') as out:
                        out.write(t.write())
                    outfiles.append(outfile)

            except AttributeError:
                pass
    for f in tolabelreg:
        for m in model_list:
            generateCtl(model=m, treefile=treefile+"."+f, seqfile=paml_msa, outfile=outfile,
                                            generateOther=False)
    t.render(treefile+".png", tree_style=ts)
    return outfiles, tolabelreg


def modelhelp():
    print ("""
    ############################################
    # MODELS
    ############################################
    modelID:    Description:             vs:
      BM     Branch model, one fg branch M0, BM with w:=0
      M0          one-ratio model;       ?
    [ M1 ]    obsolete
      M1a    branch nearly neutral       M2a
    [ M2 ]    obsolete
      M2a    branch pos. sel.            M1a
      M7     branch beta                 M8
      M8     branch beta+omega>1,alt     M7, M8a
      M8a    branch beta+omega=1,null    M8
      Ah0    branch-site                 Ah1
      Ah1    branch-site                 Ah0
    [ B ]    not yet implemented
    [ C ]    not yet implemented
    """)
    sys.exit(2)


def make_ctl_tree(treefile=None, paml_msa=None,outfile=None,
                  model="Ah0,Ah1",
                  regex=None, depth=4):
    modellist = model.split(",")
    modellist = [m.strip() for m in modellist]
    unlabeled_tree = read_tree(treefile)
    label_regex(unlabeled_tree, regex, treefile=treefile, depth=depth, model_list=modellist, paml_msa=paml_msa, outfile=outfile)


def main():

    make_ctl_tree(treefile=sys.argv[1], paml_msa=sys.argv[2], regex=sys.argv[3])

if __name__ == "__main__":
    main()