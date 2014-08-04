#!/usr/bin/python
__author__ = 'jmass'
__version__ = '0.0'

import sys
import getopt
import datetime
import os
import os.path


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
    -y, --no_copy                   only create symlink to fasta files (instead of copying)
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
    Please note that pep and nuc ids for each sequence pair needs to match exactly!
    """
    print(text)
    sys.exit(0)


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


class FastaFilesDoNotMatchException(Exception):
    pass


def check_input(dir):
    expected = [os.path.join(dir, "nuc"), os.path.join(dir, "pep")]
    d = os.path.join(dir)
    dirlist = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]
    print("checking {}".format(",".join(dirlist)))
    files_dct = {}
    tmp = None
    for e in expected:
        if e in dirlist:
            fastanames = set([".".join(name.split('.')[0:-1]) for name in os.listdir(e)])
            files_dct[e] = fastanames
            print(fastanames)
            if tmp and tmp != fastanames:
                raise FastaFilesDoNotMatchException("Fasta files in the "
                                                            "provided directories "
                                                            "do not match.\n")
            else:
                tmp = fastanames
    else:

        print("\t ok.\n")


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
        elif o in ("-y", "--no_copy"):
            no_copy = True
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
        "{:%B_%d_%Y_%H%M}".format(datetime.datetime.now())
    if not num_bootstraps:
        num_bootstraps = 100
    if not regex:
        print("No regex.\n")
        usage()

    if not labeling_level:
        labeling_level = 4
    if not models:
        models = ["Ah1", "Ah0"]
    if phase is None:
        print("No phase.\n")
        usage()
    if not num_cores:
        num_cores = 1

    check_input(input_dir)


if __name__ == "__main__":
    main()
