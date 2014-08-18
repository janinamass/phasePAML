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

@app.route('/running')
def show_running():
    return render_template('show_running.html', entries=query_db('select * from phase where status = "r" ORDER BY datetime(timestamp) DESC' ))


@app.route('/success')
def show_success():
    return render_template('show_success.html', entries=query_db('select * from phase where status ="s" ORDER BY datetime(timestamp) DESC'))


@app.route('/add', methods=['POST'])
def add_entry():
    db = get_db()
    db.execute('insert into phase (orthogroup, phase) values (?, ?)',
               [request.form['orthogroup'], request.form['phase']])
    db.commit()
    flash('New entry was successfully posted')
    return redirect(url_for('show_success'))

@app.route('/subset', methods=['POST'])
def filter_subset():
    db = get_db()
    entries = query_db('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id'] ])
    #if request.form['and']:
    #    entries = query_db('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id'] ])
    #    print('select * from phase where orthogroup = ? and run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    #elif request.form['or']:
    #    entries = query_db('select * from phase where orthogroup = ? or run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id'] ])
    #    print('select * from phase where orthogroup = ? or run_id = ? ORDER BY datetime(timestamp) DESC', [request.form['orthogroup'], request.form['run_id']])
    return show_entries(entries=entries)



if __name__ == '__main__':
    app.run(debug=True)


#todo rm all string formatting