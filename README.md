# cyc2sbml
========

Converting [Pathway Tools](http://brg.ai.sri.com/ptools/) database to Systems Biology Markup Language [SBML](http://sbml.org/).

This program reads out an organism's database from pathwaytools (e.g. metacyc) and writes a sbml model file
- Assuming that ptools is running in api mode
- Needs external python packages: [pycyc](https://github.com/ebogart/PyCyc) , [COBRApy](https://github.com/opencobra/cobrapy)


## Features
- Handling of generic metabolites
- Define substitutions (metabolites which should
- Add diffusion reactions
- Change IDs of metabolites and ractions (e.g. according to BIGG nomenclature)
-  



## Installation
First you need to install [git](http://git-scm.com/).

### [pycyc](https://github.com/ebogart/PyCyc)
```
git clone https://github.com/ebogart/PyCyc
python setup.py install
```

### [COBRApy](https://github.com/opencobra/cobrapy)
install glpk and gmp libraries (example given in debian using apt-get)
```
apt-get install libglpk-dev libgmp-dev
pip install cobra --pre
```

### cyc2sbml
```
git clone https://github.com/jotech/cyc2sbml
python cyc2sbml.py
```

## Running

Make sure that [Pathway Tools](http://brg.ai.sri.com/ptools/) is running in API mode:
```
./pathway-tools -lisp -api
```
Start cyc2sbml 
```
python cyc2sbml.py
```

### Configuration

- define reactions to be ignored in ```conf/ignore.txt```
- set diffusion reactions in ```conf/diffusion.txt```
- insert additional generic pairs in ```conf/generic_assignment.txt```
- a list of metabolites which are included even if they are generic metabolites ```conf/exceptions.txt```
- give a list of compounds which should be substituted ```conf/substitutions.txt```
