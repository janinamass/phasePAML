from flask import Flask, render_template, g, flash, request, redirect, url_for
import sqlite3
import sys
import os

DATABASE = sys.argv[1]
app = Flask(__name__)
app.secret_key = 'I am a _super_-secret key!'


def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    db.row_factory = sqlite3.Row
    return db

@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()


def query_db(query, args=(), one=False):
    cur = get_db().execute(query, args)
    rv = cur.fetchall()
    cur.close()
    return (rv[0] if rv else None) if one else rv


@app.route('/')
def home():
    cur = get_db().cursor()
    return render_template('index.html', database=os.path.abspath(DATABASE))


@app.route('/db')
def show_entries(entries=None):
    if not entries:
        return render_template('show_entries.html', entries=query_db('select * from phase ORDER BY datetime(timestamp) DESC '))
    else:
                return render_template('show_entries.html', entries=entries)
@app.route('/failed')
def show_failed():
    return render_template('show_failed.html', entries=query_db('select * from phase where status = "f" ORDER BY datetime(timestamp) DESC '))

@app.route('/running_all')
def show_running_all():
    return render_template('show_running.html', entries=query_db('select * from phase where status = "r" ORDER BY datetime(timestamp) DESC' ))

@app.route('/running')
def show_running():
    entries = query_db('select run_id, orthogroup, phase from phase '
                       'except select run_id, orthogroup, phase   from phase  where status !="r" '
                       'ORDER BY phase DESC')
    return render_template('show_running_now.html', entries=entries)


@app.route('/success')
def show_success():
    return render_template('show_success.html', entries=query_db('select * from phase where status ="s" ORDER BY datetime(timestamp) DESC'))


@app.route('/orthogroups')
def show_ortho():
    entries = query_db('select run_id, orthogroup, headers from orthoinfo')
    for e in entries:
        print(e['headers'])
        print(e['headers'].replace(',',' '))
    return render_template('show_ortho.html', entries=query_db('select * from orthoinfo'))


@app.route('/subset', methods=['POST'])
def filter_subset():
    db = get_db()
    #entries = query_db('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    #and_query = ('?', [request.form['rad_and']], )
    andor = request.form['andor']
    if andor == "AND":
        entries = query_db('select * from phase where orthogroup = ? AND  run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])

    if andor == "OR":
        entries = query_db('select * from phase where orthogroup = ? OR  run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    #or_query = ('?',[request.form['or']])
    #print(or_query)
    #if request.form['and']:
    #    entries = query_db('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id'] ])
    #    print('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    #elif request.form['or']:
    #    entries = query_db('select * from phase where orthogroup = ? or run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id'] ])
    #    print('select * from phase where orthogroup = ? or run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    return show_entries(entries=entries)

def split_space(string):
    return string.strip().split(',')


if __name__ == '__main__':
    #app.jinja_env.filters['split_space'] = split_space
    app.run(debug=True)


#select phase.*,phase.status,
# p3.status, p3.timestamp , p3.phase as p2
# from phase left join phase as p3 on phase.
# id where phase.phase=p3.phase and phase.status =
# "r" and p3.status != phase.status and phase.run_
# id=p3.run_id and phase.id > p3.id ORDER BY datetime
# (phase.timestamp) DESC;