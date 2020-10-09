__author__ = 'jmass'
import sqlite3
import sys
import subprocess
import shutil
import os
from .fastahelper import FastaParser
from .mfa2phy import mfa2phy
from .tree_labeler import make_ctl_tree
from .codeml_summary import calculatePvalue, CODEMLParser
from .map_back import main as map_back

class PipelineException(Exception):
    pass


def q(s):
    return '"' + s + '"'


def db_logger(f):
    def wrapper(*args, **kwargs):
        db = kwargs.get("db")
        run_id = kwargs.get("run_id")
        print("WRAPPER runid", run_id)
        if not run_id:
            run_id = "NULL"
        phase = kwargs.get("phase")
        status = "r"
        orthogroup = kwargs.get("orthogroup")
        connection = sqlite3.connect(db)
        with connection:
            cur = connection.cursor()
            cmd = 'INSERT INTO phase (run_id, orthogroup, phase, status) VALUES ({},{},{},{})'.format(run_id, q(orthogroup), phase,
                                                                                           q(status))
            print(cmd, db)
            cur.execute(cmd)
            connection.commit()
            try:
                res = f(*args, **kwargs)
                print(res)
                status = "s"  # success
            except PipelineException as e:
                sys.stderr.write(str(e))
                status = "f"  # fail
            cmd = 'INSERT INTO phase (run_id, orthogroup, phase, status) VALUES ({},{},{},{})'.format(run_id, q(orthogroup), phase,
                                                                                           q(status))

            print("phase {} done.\n".format(str(phase)))
            print(cmd)
            cur.execute(cmd)
            connection.commit()

    return wrapper


# phase 1: MSA pep
@db_logger
def run_prank(infile=None, outfile=None,
              cpu=1, db=None, orthogroup=None,
              run_id=None, phase=None, program=None):
    p = subprocess.Popen('{} +F -d={} -o={} '.format(program, infile, outfile),
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        #prank attaches ".best.fas"
        shutil.move(outfile + ".best.fas", outfile)
        return retval


# phase 2 pal2nal
def sort_fasta(nuc_fa=None, pep_msa=None, nuc_msa_out=None):
    if pep_msa.endswith(".msa"):
        fpp = FastaParser().read_fasta(pep_msa)
        nuc = {}
        for i in FastaParser().read_fasta(nuc_fa):
            nuc[i[0]] = i[1]
        with open(nuc_msa_out, 'w') as o:
            for f in fpp:
                o.write(">{}\n{}\n".format(f[0], nuc[f[0]]))


@db_logger
def run_pal2nal(program = None, pep_msa=None, outfile=None,
                nuc_fa=None, cpu=1, db=None,
                orthogroup=None, run_id=None,
                phase=None):
    pal2nal_nuc_in = outfile + ".nuc"
    sort_fasta(nuc_fa=nuc_fa, pep_msa=pep_msa, nuc_msa_out=pal2nal_nuc_in)
    paml = outfile + ".paml"
    pamlg = outfile + ".pamlg"
    print("CALL pal2nal nuc: {} msa: {} out:{}", pal2nal_nuc_in, pep_msa, paml)
    print("CALL pal2nal nuc: {} msa: {} out:{}", pal2nal_nuc_in, pep_msa, pamlg)

    p = subprocess.Popen('{} {} {} -output paml -nogap > {} '.format(program, pep_msa, pal2nal_nuc_in, paml),
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    #todo log errors
    print(p_out, p_err, "PAML",'{} {} {} -output paml -nogap > {} '.format(program, pep_msa, pal2nal_nuc_in, paml) )
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        p = subprocess.Popen('{} {} {} -output paml > {} '.format(program, pep_msa, pal2nal_nuc_in, pamlg),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p_out, p_err = p.communicate()
        print(p_out, p_err, "PAMLG", '{} {} {} -output paml > {} '.format(program, pep_msa, pal2nal_nuc_in, pamlg))
        retval = p.wait()
        if retval != 0:
            raise PipelineException
        elif p_err != "":
            sys.stderr.write(str(p_err))
            return retval
        else:

            return retval

@db_logger
def run_raxml(program = None, pep_msa=None, outdir=None, model=None,
              bootstrap_seed=123, num_bootstraps=None, workdir = None,
              num_cpu=None, db=None, orthogroup=None,
              run_id=None, phase=None):
    run_name = orthogroup
    pep_msa_phy = pep_msa+".phy"
    mfa2phy(pep_msa, pep_msa_phy)
    if num_bootstraps == 0:
        raxml_call= '{} -p {} -m {} -n {} -s {} -w {}'.format(program, bootstrap_seed, model, run_name, pep_msa_phy, workdir)
    else:
        raxml_call = '{} -m {} -T {} -n {} -s {} -b {} -N {} -w {}'.\
            format(program, model, num_cpu, run_name, pep_msa_phy,
                   bootstrap_seed, num_bootstraps, workdir)
    #todo raxml might have different names on other systems
    print(raxml_call)
    p = subprocess.Popen(raxml_call, shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    print(p_err)
    print(p_out)
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        if num_bootstraps == 0:
            return retval
            #todo mv to mrc (even though its not an mrc)
        else:
            bstree = os.path.join(workdir,"RAxML_bootstrap."+orthogroup)
            mrc = orthogroup+".mrc"
            # raxmlHPC -m $model -J MR -z RAxML_bootstrap.run.$nm -n $short
            raxml_call_consensus = '{} -m {} -J MR -z {} -n {} -w {}'.format(program, model, bstree, mrc, workdir)
            p = subprocess.Popen(raxml_call_consensus,
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            p_out, p_err = p.communicate()
            print(p_out)
            print(p_err)
            retval = p.wait()
            if retval != 0:
                raise PipelineException
            else:
                #rename
                return retval

@db_logger
def run_ctl_maker(paml_file=None, tree_file=None,model="Ah0,Ah1", outfile=None,
                regex=None, depth=4,
                semaphore=None, db=None,
                orthogroup=None, run_id=None,
                phase=None):
    try:
        make_ctl_tree(treefile=tree_file, paml_msa=paml_file, outfile=outfile, model=model, regex=regex, depth=depth)
    except Exception as e:
        print(e)
        raise PipelineException
    return 0
    #todo unlabeled tree from tree_file
#readtree
#labelForPamlRegex(unlabelledTree, regex, tree)
#labelForPaml(unlabelledTree,[n], tree))
#generateCtl(model=model, treefile = t, seqfile=alignment, outfile=None, generateOther = writeOther)

#needs python-qt4


@db_logger
def run_codeml(program=None, ctl_file=None, work_dir=None,
               db=None, orthogroup = None,
               run_id=None, phase=None, semaphore=None):
    codeml_call = '{} {}'.format(program, ctl_file)
    orig_wd = os.getcwd()
    os.chdir(os.path.abspath(work_dir))
    p = subprocess.Popen(codeml_call,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    print(p_out)
    print(p_err)
    retval = p.wait()
    os.chdir(orig_wd)
    if retval != 0:
        raise PipelineException
    else:
        #rename
        return retval


@db_logger
def run_seqSieve(program=None, dir=None, suffix=".msa",outdir=None,
               db=None, orthogroup=None,
               run_id=None, phase=None, semaphore=None):
    seqSieve_call = 'seqSieve -F {} -s {}'.format(dir, suffix)
    p = subprocess.Popen(seqSieve_call,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    #print(p_out)
    #print(p_err)
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        return retval


@db_logger
def run_codeml_summary(h0=None, h1=None, db=None, outfile_prefix=None, orthogroup=None,
               run_id=None, phase=None):
    try:
        pval = calculatePvalue(h0, h1)
        pval = str(pval["tablerow"])
        cp = CODEMLParser()
        sitesNEB = cp.getPositiveSitesNEB(h1)
        sitesBEB = cp.getPositiveSitesBEB(h1)
        with open(outfile_prefix+orthogroup+".csv", 'a') as out:
            out.write(pval)
        with open(outfile_prefix+"_summary.csv",'a') as out:
            out.write(pval)
        with open(outfile_prefix+orthogroup+".neb", 'a') as out:
            out.write("{}\n{}".format(str(sitesNEB["table"]), str(sitesNEB["significant"])))
        with open(outfile_prefix+orthogroup+".beb", 'a') as out:
            out.write("{}\n{}".format(str(sitesBEB["table"]), str(sitesBEB["significant"])))

    except Exception as e:
        print(e)
        raise PipelineException
    return 0


@db_logger
def run_map_back(gapped_alignment_file, list_of_positions, db=None, orthogroup=None, run_id=None, phase=None):
    try:
        retval = map_back(gappedAlignmentFile=gapped_alignment_file, listOfPositions=list_of_positions)
    except Exception:
        raise PipelineException
    if retval != 0:
        raise PipelineException
    else:
        return retval
