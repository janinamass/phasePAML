#!/usr/bin/python
__author__ = 'jmass'
__version__ = '0.1'

import sys
import getopt
import datetime
import os
import os.path
import shutil
import ConfigParser
from helpers.fastahelper import FastaParser
from helpers.dbhelper import db_check_run
from helpers.dbhelper import db_get_run_id
from helpers.wrappers import run_prank, run_pal2nal, run_raxml, run_ctl_maker, run_codeml, run_pysickle


class DirectoryExistsException(Exception):
    pass


class DifferentSequenceLengthsException(Exception):
    pass


class FastaFilesDoNotMatchException(Exception):
    pass


class HeadersDoNotMatchException(Exception):
    pass

####################################################
CONF = dict()
CONF['Directories'] = {}
CONF['Directories']['db_name'] = 'phasePAML.db'
CONF['Directories']['input_dir'] = None
CONF['Directories']['output_dir'] = None
CONF['Directories']['name'] = "{:%B_%d_%Y_%H%M}".format(datetime.datetime.now())
CONF['Pysickle'] = {}
CONF['Pysickle']['threshold'] = 9999
CONF['RAxML'] = {}
CONF['RAxML']['num_bootstraps'] = 100
CONF['RAxML']['model'] = 'PROTGAMMAJTT'
CONF['RAxML']['num_cpu'] = '8'
CONF['Labels'] = {}
CONF['Labels']['regex'] = None
CONF['Labels']['level'] = '4'
CONF['Codeml'] = {}
CONF['Codeml']['models'] = 'Ah0,Ah1'
CONF['Paths'] = {}
CONF['Paths']['prank'] = 'prank'
CONF['Paths']['raxml'] = 'raxmlHPC'
CONF['Paths']['pal2nal'] = 'pal2nal'
CONF['Paths']['codeml'] = 'codeml'
CONF['Paths']['pysickle'] = 'pysickle'
####################################################


def get_config(configfile):
    global CONF
    config = ConfigParser.ConfigParser()
    config.read(configfile)
    for section in config.sections():
        options = config.options(section)
        for option in options:
            try:
                CONF[section][option] = config.get(section, option)
            except TypeError as e:
                CONF[section] = {}
                CONF[section][option] = config.get(section, option)


def usage():
    text = """
    #############################################
    # example usage:
    # phasePAML_wrapper.py
    #        --name="testrun_42" --num_bootstraps=100  --regex='GRM*'
    #        --labeling_level=4 --models="Ah0,Ah1"
    #        --num_cores=4 --phase=0
    ############################################

    general options:
    -c, --config=CONFIG.ini         read configuration from file
    -i, --input_dir=DIR             DIR is where your input FASTA files are located,
                                    it must have a subdirectory 'nuc' and and a subdirectory 'pep'
                                    with the according files

    -o, --output_dir=DIR            DIR is where your phasePAML_run folder will be created
    -n, --name=NAME *date           identifier for the run (directory suffix, table name)
    -t, --pysickle_threshold        minimum length for paml alignment (no gaps) to pass
                                    without triggering pysickle

    -b, --num_bootstraps=INT [1000]  number of bootstraps (for RAxML runs)
    -r, --regex=REGEX               regular expression to detect nodes in the trees
                                    for foreground branch labeling

    -l, --labeling_level=INT [4]    move up INT levels from node that matches REGEX
                                    and produce a labeled tree for each of level

    -m, --models=MODELS [Ah0,Ah1]   MODELS should be a comma-separated list
                                    models: [Ah0,Ah1,M0,M1a,M2a,M7,M8,M8a,BM]

    -p, --phase=INT                 start/resume from phase INT
    -N, --num_cores=INT [1]         number of cpus to use for raxml

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


def fix_fasta(fasta=None, backupdir=None, faulty_headers=None):
    print("fasta: {}, backupdir {} ".format(fasta, backupdir))
    new_content = ""
    fp = FastaParser()
    for h, s in fp.read_fasta(fasta):
        if h in faulty_headers:
            pass
        else:
            new_content += ">{}\n{}\n".format(h, s)

    print("Info: moving original file to {}\n".format(backupdir))
    shutil.move(fasta, os.path.join(backupdir, os.path.basename(fasta)))
    with open(fasta, 'w') as fa_out:
        fa_out.write(new_content)


def check_fasta(dir, fix=False, path_dct=None):
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
    for n, p in zip(sorted(files_full_dct[expected[0]]), sorted(files_full_dct[expected[1]])):
        n = os.path.join(expected[0], n)
        p = os.path.join(expected[1], p)
        print(n, p, os.path.basename(n).split(".")[0])
        orthogroup_dct[os.path.basename(n).split(".")[0]] = []
        fpn = FastaParser().read_fasta(fasta=n)
        fpn = sorted(fpn)
        fpp = FastaParser().read_fasta(fasta=p)
        fpp = sorted(fpp)
        # check fasta length

        for i, j in zip(fpn, fpp):
            faulty =[]
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
                sys.stderr.write("Incompatible sequence length between {}:{} ({}) and {}:{} ({})".format(n, i[0], len_nuc / 3, p, j[0], len_pep))
                #raise DifferentSequenceLengthsException(
                #    "Incompatible sequence length between {}:{} ({}) and {}:{} ({})".
                #    format(n, i[0], len_nuc / 3, p, j[0], len_pep)
                #)
                faulty.append(i[0])
                #todo add pathdct
                if fix:
                    fix_fasta(n, backupdir=path_dct["faulty_nuc"], faulty_headers=faulty)
                    fix_fasta(p, backupdir=path_dct["faulty_pep"], faulty_headers=faulty)
    return orthogroup_dct


def copy_files_to_workdir(input_dir, output_dir):
    full_fasta_names = [name for name in os.listdir(input_dir)]
    for f in full_fasta_names:
        print(f, output_dir)
        shutil.copy(os.path.join(input_dir, f), os.path.join(output_dir, f))


def set_path_dct(output_dir, name):
    path_dct = {}
    base_path = os.path.join(output_dir, name)
    phase_0faulty_nuc = os.path.join(base_path, "faulty_nuc")
    phase_0faulty_pep = os.path.join(base_path, "faulty_pep")
    phase_0nuc = os.path.join(base_path, "nuc")
    phase_0pep = os.path.join(base_path, "pep")
    phase_1 = os.path.join(base_path, "MSA_pep")
    phase_2 = os.path.join(base_path, "MSA_nuc")
    phase_3 = os.path.join(base_path, "tree")
    phase_5 = os.path.join(base_path, "codeml")
    phase_6 = os.path.join(base_path, "results")
    pysickle = os.path.join(base_path, "pysickle")
    for i in [phase_0faulty_nuc, phase_0faulty_pep,
              phase_0nuc, phase_0pep, phase_1,
              phase_2, phase_3, phase_5,
              phase_6,pysickle]:
        path_dct[os.path.basename(i)] = i
    return path_dct


def make_folder_skeleton(output_dir, name):
    base_path = os.path.join(output_dir, name)
    phase_0faulty_nuc = os.path.join(base_path, "faulty_nuc")
    phase_0faulty_pep = os.path.join(base_path, "faulty_pep")
    phase_0nuc = os.path.join(base_path, "nuc")
    phase_0pep = os.path.join(base_path, "pep")
    phase_1 = os.path.join(base_path, "MSA_pep")
    phase_2 = os.path.join(base_path, "MSA_nuc")
    phase_3 = os.path.join(base_path, "tree")
    phase_5 = os.path.join(base_path, "codeml")
    phase_6 = os.path.join(base_path, "results")
    pysickle = os.path.join(base_path, "pysickle")
    if not os.path.exists(base_path):
        os.makedirs(base_path)
        for i in [phase_0faulty_nuc, phase_0faulty_pep, phase_0nuc, phase_0pep, phase_1, phase_2, phase_3, phase_5, phase_6, pysickle]:
            os.makedirs(i)
    else:
        raise DirectoryExistsException("{} already exists.".format(base_path))


def is_min_length_paml_msa(paml_msa=None, min_length=None):
    with open(paml_msa, 'r') as paml:
        line = paml.readline().strip().split(" ")
        line = [l for l in line if l.strip() != ""]
        alignment_length = int(line[1])
        print("ALIGNMENT LENGTH is {}, min length".format(alignment_length, min_length))
    if alignment_length >= min_length:
        print("true")
        return True
    else:
        print("false")
        return False


def main():
    global CONF
    configfile = None
    regex = None
    labeling_level = None
    models = None
    phase = None
    num_cores = None
    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:],
            'c:i:o:n:b:t:x:r:l:m:p:N:hHM',
            [
                'config=',
                'input_dir=',
                'output_dir=',
                'name=',
                'num_bootstraps=',
                'pysickle_threshold='
                'raxml_model='
                'regex=',
                'labeling_level=',
                'models=',
                'phase=',
                'num_cores=',
                'help',
                'HELP',
                'model_help'
            ]
        )
    except getopt.GetoptError as err:
        print (str(err))
        usage()

    for o, a in opts: #config first
        if o in ("-c", "--config"):
            configfile = a
            get_config(configfile=configfile)
    for o, a in opts:
        if o in ("-c", "--config"):
            pass
        elif o in ("-i", "--input_dir"):
            CONF['Directories']['input_dir'] = a
        elif o in ("-o", "--output_dir"):
            CONF['Directories']['output_dir'] = a
            print(CONF['Directories'])
        elif o in ("-n", "--name"):
            CONF['Directories']['name'] = a
        elif o in ("-t", "--pysickle_threshold"):
            CONF['Pysickle']['threshold'] = a
        elif o in ("-b", "--num_bootstraps"):
            CONF['RAxML']['num_bootstraps'] = a
        elif o in ("-x", "--raxml_model"):
            CONF['RAxML']['model'] = a
        elif o in ("-r", "--regex"):
            CONF['Labels']['regex'] = a
        elif o in ("-l", "--labeling_level"):
            CONF['Labels']['level'] = a
        elif o in ("-m", "--models"):
            CONF['Codeml']['models'] = a
        elif o in ("-p", "--phase"):
            phase = int(a)
        elif o in ("-N", "--num_cores"):
            CONF['RAxML']['num_cpu'] = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-H", "--HELP"):
            show_help()
        elif o in ("-M", "--model_help"):
            show_model_help()
        else:
            assert False, "unhandled option"

    try:
        input_dir = CONF['Directories']['input_dir']
    except KeyError as e:
        input_dir = None
    if not input_dir:
        print("No input directory.\n")
        usage()
    try:
        output_dir = CONF['Directories']['output_dir']
    except KeyError as e:
        output_dir = None
    if not output_dir:
        print("No output directory.\n")
        usage()
    name = CONF['Directories']['name']
    regex = CONF['Labels']['regex']
    if not regex:
        print("No regex.\n")
        usage()

    if phase is None:
        print("No phase.\n")
        usage()

    #todo if phase == 0, check source dir, otherwise check rundirs!
    #validate input folder, pep and nuc files need to match

    db = CONF["Directories"]["db_name"]
    path_dct = set_path_dct(output_dir, name=name)
    if phase == 0:
        try:
            make_folder_skeleton(output_dir=output_dir, name=name)
        except DirectoryExistsException as e:
            pass
            print(e)
            #sys.exit(1)
        try:
            orthogroup_dct = check_fasta(dir=input_dir)
        except FastaFilesDoNotMatchException as e:
            sys.stderr.write(repr(e) + '\n')
            sys.exit(1)
        except DifferentSequenceLengthsException as e:
            sys.stderr.write(repr(e) + '\n')

        copy_files_to_workdir(input_dir=os.path.join(input_dir, "nuc"), output_dir=path_dct["nuc"])
        copy_files_to_workdir(input_dir=os.path.join(input_dir, "pep"), output_dir=path_dct["pep"])
        try:
            orthogroup_dct = check_fasta(dir=os.path.join(output_dir, name), fix=True, path_dct=path_dct)
        except DifferentSequenceLengthsException as e:
            print(e)
        print(db_check_run(os.path.join(output_dir, db), run_name=name, orthogroup_dct=orthogroup_dct))
        phase = 1
    run_id = db_get_run_id(db, run_name=name)
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
                run_prank(program=CONF['Paths']['prank'], infile=infile,
                          outfile=outfile,
                          db=db,
                          run_id=run_id,
                          orthogroup=orthogroup,
                          phase=1)
        phase = 2
    if phase == 2:
        ###phase 2:
        nucfiles = os.listdir(path_dct["nuc"])
        nuc_msa_dct = {}
        for f in os.listdir(path_dct["MSA_pep"]):
            if f.endswith(".msa"): #ignore .reduced
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
            #if not pep_msa.endswith('.msa'):
            #    continue
            #produces .pamlg, .paml
            run_pal2nal(program=CONF['Paths']['pal2nal'], pep_msa=pep_msa, nuc_fa=nuc_fa,
                        outfile=os.path.join(path_dct["MSA_nuc"],
                                             orthogroup),
                        cpu=1,
                        db=db,
                        orthogroup=orthogroup,
                        run_id=run_id,
                        phase=2)
        # copy .paml files to codeml dir
        for pamlfile in os.listdir(path_dct["MSA_nuc"]):
            if pamlfile.endswith(".paml"):
                shutil.copy(os.path.join(path_dct["MSA_nuc"], pamlfile), os.path.join(path_dct["codeml"], pamlfile))
                #check if the files are min length, copy pep msa file to pysickle if they're not
                if not is_min_length_paml_msa(os.path.join(path_dct["MSA_nuc"], pamlfile), min_length=int(CONF['Pysickle']['threshold'])):
                    #todo fix qnd hack
                    orthogroup = pamlfile.split(".")[0]
                    shutil.copy(os.path.join(path_dct["MSA_pep"], orthogroup+".msa"), os.path.join(path_dct["pysickle"], orthogroup+".msa"))
                    nuc = [n for n in os.listdir(path_dct["nuc"]) if n.split(".")[0] == orthogroup][0]
                    shutil.copy(os.path.join(path_dct["nuc"], nuc), os.path.join(path_dct["pysickle"], nuc))
                    #todo if pysickle
        pysickle=True
        if pysickle:
            print("running pysickle")
            run_pysickle(program=CONF['Paths']['pysickle'], dir=path_dct["pysickle"], db=db,orthogroup="__pysickle__",run_id=run_id,phase=999)
            for pysickled in os.listdir(os.path.join(path_dct["pysickle"], "ps_out_si")):
                if pysickled.endswith(".tmp"):#marks new pep.fa
                    print(pysickled)
                    new_name = pysickled.replace(".", "_").replace("_tmp", ".fa")
                    print(new_name, "new name", "infile:", os.path.join(path_dct["pysickle"],"ps_out_si", new_name))
                    run_prank(program=CONF['Paths']['prank'],infile=os.path.join(path_dct["pysickle"],"ps_out_si", pysickled),
                              outfile= os.path.join(path_dct["MSA_pep"], new_name.split(".")[0]+".msa"),
                              cpu=1, db=db,
                              orthogroup=new_name.split(".")[0], run_id=run_id,
                              phase=99)
                    #nwo we still need the new nuc
                    nucfa = [n for n in os.listdir(path_dct["nuc"]) if n.split(".")[0] == pysickled.split(".")[0]][0]
                    orthogroup=new_name.split(".")[0]
                    outfile = os.path.join(path_dct["MSA_nuc"], orthogroup+".msa")
                    run_pal2nal(program=CONF['Paths']['pal2nal'], pep_msa=os.path.join(path_dct["MSA_pep"], new_name.split(".")[0]+".msa"),
                                outfile=outfile, nuc_fa=os.path.join(path_dct["MSA_nuc"],nucfa), db=db, phase=10,run_id=run_id, orthogroup=orthogroup)
                    #merge back
                    shutil.copy(outfile+".paml", os.path.join(path_dct["codeml"], orthogroup+".paml"))

        phase = 3 #raxml
    if phase == 3: #this will fail if raxml was run there before
        for pep_msa in os.listdir(path_dct["MSA_pep"]):
            if pep_msa.endswith(".msa"):
                orthogroup = os.path.basename(pep_msa).split(".")[0]
                run_raxml(program=CONF['Paths']['raxml'], pep_msa=os.path.join(path_dct["MSA_pep"], pep_msa), outdir=path_dct['tree'],
                          num_bootstraps=int(CONF['RAxML']['num_bootstraps']), db=db, model=CONF['RAxML']['model'],
                          orthogroup=orthogroup, run_id=run_id,
                          phase=phase, workdir=os.path.abspath(path_dct["tree"]))
        phase = 4

    if phase == 4:
        for mrc in os.listdir(path_dct["tree"]):
            if mrc.endswith(".mrc") and mrc.startswith("RAxML_MajorityRuleConsensusTree."):
                new_name = mrc.split("RAxML_MajorityRuleConsensusTree.")[1]
                shutil.copy(os.path.join(path_dct["tree"], mrc), os.path.join(path_dct["codeml"], new_name))
            elif mrc.startswith("RAxML_result"): #non bootstrap trees
                new_name = mrc.split("RAxML_result.")[1]+".mrc"
                shutil.copy(os.path.join(path_dct["tree"], mrc), os.path.join(path_dct["codeml"], new_name))
        #todo remove tree_lab folder
        for mrc in os.listdir(path_dct["codeml"]):
            if mrc.endswith(".mrc"):
                orthogroup = mrc.split(".")[0]
                pamlfile = os.path.join(path_dct["codeml"], orthogroup+".paml")
                treefile = os.path.join(path_dct["codeml"], orthogroup+".mrc")
                print(orthogroup)
                print(pamlfile, mrc)
                print(regex)

                run_ctl_maker(paml_file=pamlfile, tree_file=treefile, model=CONF['Codeml']['models'],
                              outfile=treefile, regex=CONF['Labels']['regex'],
                              depth=int(CONF['Labels']['level']), db=db,
                              orthogroup=orthogroup, run_id=run_id,
                              phase=phase)
        phase = 5
    if phase == 5:
        for ctl in os.listdir(path_dct["codeml"]):
            if ctl.endswith(".ctl"):
                workdir = path_dct["codeml"]
                orthogroup = ctl.split(".")[0]
                run_codeml(program=CONF['Paths']['codeml'], ctl_file=ctl, work_dir=workdir,
                           db=db, orthogroup=orthogroup,
                           run_id=run_id, phase=phase)


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

#todo if paml < X, cp to pysickle
#todo add pysickle to getopt: on length <x, always, never, default:never