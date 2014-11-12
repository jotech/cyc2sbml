# cyc2sbml

Converting [Pathway Tools](http://brg.ai.sri.com/ptools/) database to Systems Biology Markup Language [SBML](http://sbml.org/).

This program reads out an organism's database from pathwaytools (e.g. metacyc) and writes a sbml model file.
- Assuming that Pathway Tools is running in api mode
- Needs external python packages: [pycyc](https://github.com/ebogart/PyCyc) , [cobrapy](https://github.com/opencobra/cobrapy)


## Features
- Automatical handling of generic metabolites (these are unspecific or abstract compounds which have otherweise to be specified manually. For example NADH-P-OR-NOP in Pathway Tools which could stand for NADH and NADPH)
- Define substitutions (metabolites which should be different in sbml)
- Add standard diffusion reactions (water, ethanol and small, uncharged <=C3 metabolites can pass the membrane)
- Change IDs and names of metabolites and reactions (e.g. according to the BIGG nomenclature)
- Ignore reaction that are not important in the metabolism)
- Change Pathway Tools gene names


## Installation
First you need to install [git](http://git-scm.com/).

### [pycyc](https://github.com/ebogart/PyCyc)
```
git clone https://github.com/ebogart/PyCyc
python setup.py install
```

### [cobrapy](https://github.com/opencobra/cobrapy)
install glpk and gmp libraries (example given in debian using apt-get)
```
apt-get install libglpk-dev libgmp-dev
pip install cobra --pre
```

### cyc2sbml
```
git clone https://github.com/jotech/cyc2sbml
```

## Excecution

Make sure that [Pathway Tools](http://brg.ai.sri.com/ptools/) is running in API mode:
```
./pathway-tools -lisp -api
```
Start cyc2sbml (no parameters needed, all important decisions will clarified interactively)
```
python cyc2sbml.py
```

### Configuration
Settings for all configurations are given and can be switched on/off interactively during runtime.

- define reactions to be ignored in ```conf/ignore.txt```
- set diffusion reactions in ```conf/diffusion.txt```
- insert additional generic pairs in ```conf/generic_assignment.txt```
- a list of metabolites which are included even if they are generic metabolites ```conf/exceptions.txt```
- give a list of compounds which should be substituted ```conf/substitutions.txt```
