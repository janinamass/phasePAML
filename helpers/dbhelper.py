import sqlite3
import os


def q(s):
    return '"' + s + '"'


def check_phase(cursor, run_name):
    pass


def db_check_run(db, run_name, orthogroup_dct):
    """check if name of run already exists in table"""
    if os.path.isfile(db):
        print("{} exists".format(db))
        con = sqlite3.connect(db)
        cmd = "SELECT * FROM run WHERE name = " + q(run_name)
        with con:
            cur = con.cursor()
            cur.execute(cmd)
            res = cur.fetchall()
        if len(res) > 0:
            print("Run {} already exists.\n".format(run_name))
            check_phase(run_name)
    else:
        print("Creating new db file {}.\n".format(db))
        cmd = 'CREATE TABLE run (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT);'
        con = sqlite3.connect(db)
        with con:
            cur = con.cursor()
            cur.execute(cmd)
            con.commit()
            cmd = "INSERT INTO run(name) VALUES({});".format(q(run_name))
            cur.execute(cmd)
            con.commit()
            cmd = "SELECT id FROM run WHERE name = " + q(run_name)
            cur.execute(cmd)
            run_id = cur.fetchone()[0]
            init_orthoinfotable(connection=con,
                                run_id=run_id,
                                orthogroup_dct=orthogroup_dct
            )

            init_phasetable(connection=con,
                            run_id=run_id,
                            orthogroup_list=orthogroup_dct.keys()
            )

            con.commit()


def init_phasetable(connection, run_id, orthogroup_list):
    phase = 0
    cmd = 'CREATE TABLE phase (' \
          'id INTEGER PRIMARY KEY AUTOINCREMENT, ' \
          'run_id INTEGER, ' \
          'phase INTEGER, ' \
          'orthogroup TEXT, ' \
          'FOREIGN KEY(run_id) REFERENCES run(id), ' \
          'FOREIGN KEY(orthogroup) REFERENCES orthoinfo(orthogroup)' \
          ');'
    cur = connection.cursor()
    cur.execute(cmd)
    print("phase table created.\n")
    cmdf = lambda runid, orthgroup: \
        'INSERT INTO phase(run_id, phase, orthogroup)' \
        ' VALUES ({},{},{})'.format(runid, phase, orthgroup)

    for o in orthogroup_list:
        o = q(o)

        cmd = cmdf(run_id, o)
        print(cmd, "CMD")
        print(cmd)
        cur.execute(cmd)


def init_orthoinfotable(connection, run_id, orthogroup_dct):
    cmd = 'CREATE TABLE orthoinfo (' \
          'run_id INTEGER, ' \
          'orthogroup TEXT, ' \
          'headers TEXT,' \
          'FOREIGN KEY(run_id) REFERENCES run(id), ' \
          'PRIMARY KEY(run_id, orthogroup)' \
          ');'
    cur = connection.cursor()
    cur.execute(cmd)
    print("phase table created.\n")
    connection.commit()
    cmdf = lambda runid, orthgroup, headers: \
        'INSERT INTO orthoinfo(run_id, orthogroup, headers)' \
        ' VALUES ({},{},{})'.format(runid, orthgroup, headers)

    for o in orthogroup_dct:
        qo = q(o)  # orthogroup_name
        cmd = cmdf(run_id, qo, q(",".join(orthogroup_dct[o])))
        print(cmd)
        cur.execute(cmd)