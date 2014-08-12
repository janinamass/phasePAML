#!/usr/bin/python
__author__ = 'jmass'
__version__ = '0.0'

import sys
import getopt
import datetime
import os
import os.path
import shutil

from helpers.fastahelper import FastaParser
from helpers.dbhelper import db_check_run
from helpers.dbhelper import db_get_run_id
from helpers.wrappers import run_prank, run_pal2nal, run_raxml, run_ctl_maker


class DirectoryExistsException(Exception):
    pass


class DifferentSequenceLengthsException(Exception):
    pass


class FastaFilesDoNotMatchException(Exception):
    pass


class HeadersDoNotMatchException(Exception):
    pass



def usage():
    text = """
    #############################################
    # example usage:
    # phasePAML_wrapper.py
    #        --name="testrun_42" --num_bootstraps=100  --regex='GRM.*'
    #        --labeling_level=4 --models="Ah0,Ah1"
    #        --num_cores=4 --phase=0 --no_copy
    ############################################

    general options:
    -i, --input_dir=DIR             DIR is where your input FASTA files are located,
                                    it must have a subdirectory 'nuc' and and a subdirectory 'pep'
                                    with the according files

    -o, --output_dir=DIR            DIR is where your phasePAML_run folder will be created
    -n, --name=NAME *date           identifier for the run (directory suffix, table name)
    -b, --num_bootstraps=INT [100]  number of bootstraps (for RAxML runs)
    -r, --regex=REGEX               regular expression to detect nodes in the trees
                                    for foreground branch labeling

    -l, --labeling_level=INT [4]    move up INT levels from node that matches REGEX
                                    and produce a labeled tree for each of level

    -m, --models=MODELS [Ah0,Ah1]   MODELS should be a comma-separated list
                                    models: [Ah0,Ah1,M0,M1a,M2a,M7,M8,M8a,BM]

    -p, --phase=INT                 start/resume from phase INT
    -c, --num_cores=INT [1]         number of cpus to use
    # -y, --no_copy                   only create symlink to fasta files (instead of copying)
    -h, --help                      prints this
    -H, --HELP                      more help
    -M, --model_help                more help on models
    """

    print(text)
    sys.exit(2)


def show_help():
    # TODO add extended help message
    text = """
    #############################################
    #
    #
    #
    #
    #############################################
    Please note that pep and nuc ids for each sequence pair need to match exactly!
    """
    print(text)
    sys.exit(0)


def show_phase_help():
    # TODO detailed phase description with dependencies
    text = """
    ##################################################################
    # phase         description      dependencies
    # 0             from scratch     valid pep and nuc files
    # 1             MSA              prank
    # 2             nucMSA           pal2nal
    # 3             tree             RAxML
    # 4             tree labeling    ete2
    # 5             select.analysis  codeml (PAML)
    # 6             summarize        result from phase 4 for all files
    ###################################################################
    """


def show_model_help():
    text = """
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
    """
    print(text)
    sys.exit(0)


def check_fasta(dir):
    expected = [os.path.join(dir, "nuc"), os.path.join(dir, "pep")]
    d = os.path.join(dir)
    dirlist = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]
    tmp = None
    files_dct = {}
    files_full_dct = {}
    for e in expected:
        if e in dirlist:
            full_fasta_names = [name for name in os.listdir(e)]
            fasta_names = set([".".join(name.split('.')[0:-1]) for name in os.listdir(e)])
            files_full_dct[e] = full_fasta_names
            if tmp and tmp != fasta_names:
                raise FastaFilesDoNotMatchException("Fasta files in the "
                                                    "provided directories "
                                                    "do not match.")
            else:
                tmp = fasta_names
    orthogroup_dct = {}
    for n, p in zip(files_full_dct[expected[0]], files_full_dct[expected[1]]):
        n = os.path.join(expected[0], n)
        p = os.path.join(expected[1], p)
        print(n, p, os.path.basename(n).split(".")[0] )
        orthogroup_dct[os.path.basename(n).split(".")[0]] = []
        fpn = FastaParser().read_fasta(fasta=n)
        fpn = sorted(fpn)
        fpp = FastaParser().read_fasta(fasta=p)
        fpp = sorted(fpp)
        # check fasta length
        for i, j in zip(fpn, fpp):
            orthogroup_dct[os.path.basename(n).split(".")[0]].append(i[0])
            #print(i[0], j[0])
            if i[0] != j[0]:
                raise HeadersDoNotMatchException(
                    "Headers do not match between your pep and nuc files.\n"
                    "{}:{} vs {}:{}".format(n, i[0], p, j[0]))
            len_nuc = len(i[1])
            len_pep = len(j[1])
            print(len_nuc, len_nuc / 3, len_pep)
            if not len_pep == len_nuc / 3 \
                    and not (len_pep - 1 == len_nuc / 3 and j[1][-1] == "*") \
                    and not (len_pep + 1 == len_nuc / 3 and j[1][-1] != "*"):
                raise DifferentSequenceLengthsException(
                    "Incompatible sequence length between {}:{} ({}) and {}:{} ({})".
                    format(n, i[0], len_nuc / 3, p, j[0], len_pep)
                )
    print(orthogroup_dct)
    return orthogroup_dct


def copy_files_to_workdir(input_dir, output_dir):
    full_fasta_names = [name for name in os.listdir(input_dir)]
    for f in full_fasta_names:
        print(f, output_dir)
        shutil.copy(os.path.join(input_dir, f), os.path.join(output_dir, f))


def set_path_dct(output_dir, name):
    path_dct = {}
    base_path = os.path.join(output_dir, name)
    phase_0nuc = os.path.join(base_path, "nuc")
    phase_0pep = os.path.join(base_path, "pep")
    phase_1 = os.path.join(base_path, "MSA_pep")
    phase_2 = os.path.join(base_path, "MSA_nuc")
    phase_3 = os.path.join(base_path, "tree")
    phase_4 = os.path.join(base_path, "tree_lab")
    phase_5 = os.path.join(base_path, "codeml")
    phase_6 = os.path.join(base_path, "results")
    for i in [phase_0nuc, phase_0pep, phase_1,
              phase_2, phase_3, phase_4, phase_5,
              phase_6]:
        path_dct[os.path.basename(i)] = i
    return path_dct


def make_folder_skeleton(output_dir, name):
    base_path = os.path.join(output_dir, name)
    phase_0nuc = os.path.join(base_path, "nuc")
    phase_0pep = os.path.join(base_path, "pep")
    phase_1 = os.path.join(base_path, "MSA_pep")
    phase_2 = os.path.join(base_path, "MSA_nuc")
    phase_3 = os.path.join(base_path, "tree")
    phase_4 = os.path.join(base_path, "tree_lab")
    phase_5 = os.path.join(base_path, "codeml")
    phase_6 = os.path.join(base_path, "results")
    if not os.path.exists(base_path):
        os.makedirs(base_path)
        for i in [phase_0nuc, phase_0pep, phase_1, phase_2, phase_3, phase_4, phase_5, phase_6]:
            os.makedirs(i)
    else:
        raise DirectoryExistsException("{} already exists.".format(base_path))

DB = "phasePAML.db"


def main():
    input_dir = None
    output_dir = None
    name = None
    num_bootstraps = None
    regex = None
    labeling_level = None
    models = None
    phase = None
    num_cores = None
    no_copy = False
    #todo rm no_copy option
    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:],
            'i:o:n:b:r:l:m:p:c:yhHM',
            [
                'input_dir=',
                'output_dir=',
                'name=',
                'num_bootstraps=',
                'regex=',
                'labeling_level=',
                'models=',
                'phase=',
                'num_cores=',
                'no_copy',
                'help',
                'HELP',
                'model_help'
            ]
        )
    except getopt.GetoptError as err:
        print (str(err))
        usage()

    for o, a in opts:
        if o in ("-i", "--input_dir"):
            input_dir = a
        elif o in ("-o", "--output_dir"):
            output_dir = a
        elif o in ("-n", "--name"):
            name = a
        elif o in ("-b", "--num_bootstraps"):
            num_bootstraps = int(a)
        elif o in ("-r", "--regex"):
            regex = a
        elif o in ("-l", "--labeling_level"):
            labeling_level = int(a)
        elif o in ("-m", "--models"):
            models = a
        elif o in ("-p", "--phase"):
            phase = int(a)
        elif o in ("-c", "--num_cores"):
            num_cores = int(a)
        #elif o in ("-y", "--no_copy"):
            #no_copy = True
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-H", "--HELP"):
            show_help()
        elif o in ("-M", "--model_help"):
            show_model_help()
        else:
            assert False, "unhandled option"

    if not input_dir:
        print("No input directory.\n")
        usage()
    if not output_dir:
        print("No output directory.\n")
        usage()
    if not name:
        name = "{:%B_%d_%Y_%H%M}".format(datetime.datetime.now())
    if not num_bootstraps:
        num_bootstraps = 100
    if not regex:
        print("No regex.\n")
        usage()

    if not labeling_level:
        labeling_level = 4
    if not models:
        models = "Ah1,Ah0"
    if phase is None:
        print("No phase.\n")
        usage()
    if not num_cores:
        num_cores = 1
    # ##################
    # todo sanity check:
    #prank
    #mafft
    #raxml
    #pal2nal
    #codeml
    ################### raise Warning, but offer to continue [y,N]V

    #todo if phase == 0, check source dir, otherwise check rundirs!
    #validate input folder, pep and nuc files need to match
    path_dct = set_path_dct(output_dir, name=name)
    if phase == 0:
        try:
            orthogroup_dct = check_fasta(dir=input_dir)
        except (FastaFilesDoNotMatchException, DifferentSequenceLengthsException) as e:
            sys.stderr.write(repr(e) + '\n')
            sys.exit(1)
        try:
            make_folder_skeleton(output_dir=output_dir, name=name)
        except DirectoryExistsException as e:
            pass
            print(e)
            #sys.exit(1)
        copy_files_to_workdir(input_dir=os.path.join(input_dir, "nuc"), output_dir=path_dct["nuc"])
        copy_files_to_workdir(input_dir=os.path.join(input_dir, "pep"), output_dir=path_dct["pep"])
        print(orthogroup_dct, "OGD")
        print(db_check_run(os.path.join(output_dir, DB), run_name=name, orthogroup_dct=orthogroup_dct))
        phase = 1
    run_id = db_get_run_id(DB, run_name=name)
    if not run_id:
        sys.stderr.write("Got No run_id, check your db.\n")
        #sys.exit(1)
    elif len(run_id) > 1:
        sys.stderr.write("Got multiple run_ids with the same name, check your db.\n")
        #sys.exit(1)
    else:
        run_id = run_id[0][0]
    print("Info: run_id is {}".format(run_id))

    if phase == 1:
        if not os.path.exists(path_dct["pep"]):
            sys.stderr.write("Could not find {}".format(path_dct["pep"]))
            sys.exit(1)
        else:
            for f in os.listdir(path_dct["pep"]):
                infile = os.path.join(path_dct["pep"], f)
                outfile = os.path.join(path_dct["MSA_pep"], f.split(".")[0]+".msa")
                orthogroup = os.path.basename(infile).split(".")[0]
                run_prank(infile=infile,
                          outfile=outfile,
                          db=DB,
                          run_id=run_id,
                          orthogroup=orthogroup,
                          phase=1)
        phase = 2
    if phase == 2:
        ###phase 2:
        nucfiles = os.listdir(path_dct["nuc"])
        nuc_msa_dct = {}
        for f in os.listdir(path_dct["MSA_pep"]):
            msa_file = os.path.join(path_dct["MSA_pep"], f)
            nucfile = [os.path.join(path_dct["nuc"], n) for n in nucfiles if n.split(".")[0] == f.split(".")[0]]
            if len(nucfile) > 1:
                sys.stderr.write("Sth wrong with your fasta files,"
                                " got multiple matches for {}".format(f.split(".")[0]))
                sys.exit(1)
            else:
                nucfile = nucfile[0]
                nuc_msa_dct[nucfile] = msa_file
        #only pass if all files have a match
        for nuc_fa, pep_msa in nuc_msa_dct.items():
            orthogroup = os.path.basename(pep_msa).split(".")[0]
            #produces .pamlg, .paml
            run_pal2nal(pep_msa=pep_msa, nuc_fa=nuc_fa,
                        outfile=os.path.join(path_dct["MSA_nuc"],
                                             orthogroup),
                        cpu=1,
                        semaphore=None,
                        db=DB,
                        orthogroup=orthogroup,
                        run_id=run_id,
                        phase=2)
        # copy .paml files to codeml dir
        for pamlfile in os.listdir(path_dct["MSA_nuc"]):
            if pamlfile.endswith(".paml"):
                shutil.copy(os.path.join(path_dct["MSA_nuc"], pamlfile), os.path.join(path_dct["codeml"], pamlfile))
        phase = 3 #raxml
    if phase == 3:
        for pep_msa in os.listdir(path_dct["MSA_pep"]):
            if pep_msa.endswith(".msa"):
                orthogroup = os.path.basename(pep_msa).split(".")[0]
                #todo raxml model, for now always PROTGAMMAJTT
                run_raxml(pep_msa=os.path.join(path_dct["MSA_pep"], pep_msa), outdir=path_dct['tree'],
                          num_bootstraps=5, db=DB, model="PROTGAMMAJTT",
                          orthogroup=orthogroup, run_id=run_id,
                          phase=phase, workdir=os.path.abspath(path_dct["tree"]))
        phase = 4

    if phase == 4:
        #copy majority rule consensus trees to tree_lab folder
        for mrc in os.listdir(path_dct["tree"]):
            if mrc.endswith(".mrc") and mrc.startswith("RAxML_MajorityRuleConsensusTree."):
                new_name = mrc.split("RAxML_MajorityRuleConsensusTree.")[1]
                shutil.copy(os.path.join(path_dct["tree"], mrc), os.path.join(path_dct["codeml"], new_name))
        #todo remove tree_lab folder
        for mrc in os.listdir(path_dct["codeml"]):
            if mrc.endswith(".mrc"):
                orthogroup = mrc.split(".")[0]
                pamlfile = os.path.join(path_dct["codeml"], orthogroup+".paml")
                treefile =  os.path.join(path_dct["codeml"], orthogroup+".mrc")
                print(orthogroup)
                print(pamlfile, mrc)
                print(regex)
                run_ctl_maker(paml_file=pamlfile, tree_file=treefile, model=models,
                              outfile=treefile, regex=regex,
                              depth=4, semaphore=None, db=DB,
                              orthogroup=orthogroup, run_id=run_id,
                              phase=phase)



    # raise OverwriteRunException if started with phase below completed phases
    #catch OverwriteRunException :"{} has completed phase n, but you want to start phase m<n.
    # Are you sure you want to overwrite phases x,yz [y,N]?"
    # if overwrite: drop rows from table
    #


if __name__ == "__main__":
    main()
    #todo add run to database if db already exists
    #todo orthogroups in table should not have file ext
    #todo kick out corrupted fasta
    #todo enable resume

