import sqlite3
import os


def q(s):
    return '"' + s + '"'


def check_phase(cursor=None, run_name=None):
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
                check_phase(cur, run_name)
            else:
                cmd = "INSERT INTO run(name) VALUES({});".format(q(run_name))
                cur.execute(cmd)
                con.commit()
                cmd = "SELECT id FROM run WHERE name = " + q(run_name)
                cur.execute(cmd)
                run_id = cur.fetchone()[0]
                init_orthoinfotable(connection=con,
                                    run_id=run_id,
                                    orthogroup_dct=orthogroup_dct)

                init_phasetable(connection=con,
                                run_id=run_id,
                                orthogroup_list=orthogroup_dct.keys())

                con.commit()


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
                                orthogroup_dct=orthogroup_dct)

            init_phasetable(connection=con,
                            run_id=run_id,
                            orthogroup_list=orthogroup_dct.keys())

            con.commit()


def init_phasetable(connection, run_id, orthogroup_list):
    cmd = 'CREATE TABLE phase (' \
          'id INTEGER PRIMARY KEY AUTOINCREMENT, ' \
          'run_id INTEGER, ' \
          'phase INTEGER, ' \
          'orthogroup TEXT, ' \
          'status TEXT, ' \
          'timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,' \
          'FOREIGN KEY(run_id) REFERENCES run(id), ' \
          'FOREIGN KEY(orthogroup) REFERENCES orthoinfo(orthogroup)' \
          ');'
    cur = connection.cursor()
    cur.execute(cmd)
    print("phase table created.\n")
    # status
    # f: failed
    # r: running
    # s: success
    cmdf = lambda runid, phase, orthogroup, status: \
        'INSERT INTO phase(run_id, phase, orthogroup, status)' \
        ' VALUES ({},{},{},{})'.format(runid, phase, orthogroup, status)

    phase = 0
    for o in orthogroup_list:
        cmd = cmdf(run_id, phase, q(o), q("s"))
        print(cmd, "CMD")
        print(cmd)
        cur.execute(cmd)
    connection.commit()


def update_phasetable(db, run_id, orthogroup_id, phase, status):
    # todo test
    # update might be misleading, just insert another line for the current step
    # status
    # f: failed
    # r: running
    # s: success
    con = sqlite3.connect(db)
    with con:
        cur = con.cursor()
        cmdf = lambda runid, phase, orthogroup_id, status: \
            'INSERT INTO phase(run_id, phase, orthogroup, status)' \
            ' VALUES ({},{},{},{})'.format(runid, phase, orthogroup_id, status)

        cmd = cmdf(run_id, phase, q(orthogroup_id), q(status))
        print(cmd, "CMD")
        cur.execute(cmd)
        con.commit()


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
    print("orthoinfo table created.\n")
    connection.commit()
    cmdf = lambda runid, orthgroup, headers: \
        'INSERT INTO orthoinfo(run_id, orthogroup, headers)' \
        ' VALUES ({},{},{})'.format(runid, orthgroup, headers)

    for o in orthogroup_dct:
        print(o, "OO")
        qo = q(o)  # orthogroup_n   ame
        cmd = cmdf(run_id, qo, q(",".join(orthogroup_dct[o])))
        print(cmd)
        cur.execute(cmd)
    connection.commit()


def db_get_run_id(db, run_name):
    con = sqlite3.connect(db)
    with con:
        cur = con.cursor()
        cmd = "SELECT id FROM run WHERE name = {};".format(q(run_name))
        print(cmd)
        cur.execute(cmd)
        run_ids = cur.fetchall()
    return run_ids