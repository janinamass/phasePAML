from flask import Flask, render_template
import sqlite3
from flask import g
from flask import render_template
import sys

DATABASE = sys.argv[1]
#todo set db via cmd line
app = Flask(__name__)


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
    return render_template('index.html')

@app.route('/db')
def show_entries():
    db = get_db()
    res = []
    for orthogroup in query_db('select * from phase'):
        res.append("{} {} {}".format(orthogroup['orthogroup'], 'has the id', orthogroup['id']))
    print(res)
    return render_template('show_entries.html', entries=res)

if __name__ == '__main__':
    app.run(debug=True)
    for orthogroup in query_db('select * from phase'):
        print orthogroup['orthogroup'], 'has the id', orthogroup['id']


#todo rm all string formatting