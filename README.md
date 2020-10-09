phasePAML
=========




additionally needs (python): 

    # python3 -m venv env
    # source env/bin/activate
    (env) pip install ete3 
    (env) pip install matplotlib
    (env) pip install numpy
    (env) pip install seqsieve 
    (env) pip install PyQt5
    (env) pip install scipy
external:

     - pal2nal
     - prank
     - raxml
     - paml
    
run with:

    python phasePAML_wrapper.py -c myconfig.ini -i myinput/ -o myoutput -p0
    
    
