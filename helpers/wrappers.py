__author__ = 'jmass'
import sqlite3
import sys
import subprocess
import shutil
from fastahelper import FastaParser

class PipelineException(Exception):
    pass


def q(s):
    return '"' + s + '"'


def db_logger(f):
    def wrapper(*args, **kwargs):
        db = kwargs.get("db")
        run_id = kwargs.get("run_id")
        phase = kwargs.get("phase")
        status = "r"
        orthogroup = kwargs.get("orthogroup")
        connection = sqlite3.connect(db)
        with connection:
            cur = connection.cursor()
            cmd = 'INSERT INTO phase (orthogroup, phase, status) VALUES ({},{},{})'.format(q(orthogroup), phase, q(status))
            print(cmd)
            cur.execute(cmd)
            connection.commit()
            try:
                res = f(*args, **kwargs)
                print(res)
                status = "s" # success
            except PipelineException as e:
                sys.stderr.write(str(e))
                status = "f" # fail
            cmd = 'INSERT INTO phase (orthogroup, phase, status) VALUES ({},{},{})'.format(q(orthogroup), phase, q(status))

            print("phase {} done.\n".format(str(phase)))
            print(cmd)
            cur.execute(cmd)
            connection.commit()
    return wrapper


# phase 1: MSA pep
@db_logger
def run_prank(infile=None, outfile=None,
              cpu=1, semaphore=None,
              db=None, orthogroup=None,
              run_id=None, phase=None):
    p = subprocess.Popen('prank +F -d={} -o={}'.format(infile, outfile),
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    #p_out, p_err = p.communicate()
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        #prank attaches ".best.fas"
        shutil.move(outfile+".best.fas", outfile)
        return retval


#phase 2 pal2nal
def sort_fasta(nuc_fa=None, pep_msa=None, nuc_msa_out=None):
    fpp = FastaParser().read_fasta(pep_msa)
    nuc = {}
    for i in FastaParser().read_fasta(nuc_fa):
        nuc[i[0]] = i[1]
    with open(nuc_msa_out, 'w') as o:
        for f in fpp:
            o.write(">{}\n{}\n".format(f[0], nuc[f[0]]))


@db_logger
def run_pal2nal(pep_msa=None, outfile=None,
                nuc_fa=None, cpu=1,
                semaphore=None,
                db=None, orthogroup=None,
                run_id=None, phase=None):
    pal2nal_nuc_in = outfile+".nuc"
    sort_fasta(nuc_fa=nuc_fa, pep_msa=pep_msa, nuc_msa_out=pal2nal_nuc_in)
    paml=outfile+".paml"
    pamlg = outfile+".pamlg"
    print("CALL pal2nal nuc: {} msa: {} out:{}", pal2nal_nuc_in, pep_msa, paml)
    print("CALL pal2nal nuc: {} msa: {} out:{}", pal2nal_nuc_in, pep_msa, pamlg)

    p = subprocess.Popen('pal2nal.pl {} {} -output paml -nogap > {} '.format(pep_msa, nuc_fa, paml),
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p_out, p_err = p.communicate()
    #print(p_out, p_err)
    retval = p.wait()
    if retval != 0:
        raise PipelineException
    else:
        p = subprocess.Popen('pal2nal.pl {} {} -output paml > {} '.format(pep_msa, nuc_fa, pamlg),
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
        p_out, p_err = p.communicate()
        #print(p_out, p_err)
        retval = p.wait()
        if retval != 0:
            raise PipelineException
        else:

            return retval



@db_logger
def run_raxml(infile=None, outfile=None,
              cpu=1, semaphore=None,
              db=None, orthogroup=None,
              run_id=None, phase=None,
              status=None):
    i = 0
    for i in range(0, 100):
        print(i)
        i += 1

    return i
#todo nr bootstraps etc
