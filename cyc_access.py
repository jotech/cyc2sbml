import re
from cobra import Model, Reaction, Metabolite


def test():
  print "hello"


def reaction_name(org, reaction):
  """Returns the full name of a reaction"""
  return org.get_name_string(reaction)


def reaction_subsystem(org, reaction):
  """returns subsystem of a reaction"""
  pathways = reaction_pathways(org, reaction)
  if len(pathways) > 0:
    return str(pathways[0])
  else:
    return str(org.reaction_type(reaction))
   


def no_style(string):
  """Get rid of annoying metacyc style in strings"""
  string = re.sub("/","",string)
  string = string.replace("<sub>", "").replace("<SUB>", "").replace("<sup>", "").replace("<SUP>", "")
  string = string.replace("<i>", "").replace("<SUB>", "").replace("<sup>", "").replace("<SUP>", "")
  return string


def reaction_pathways(org, reaction):
  """Returns a list of all pathways""" 
  pathways = []
  pws = reaction.in_pathway
  if pws != None:
    if hasattr(pws, '__iter__'): # check if iterable
      i = 0
      for pw in pws:
        if str(pw).find("RXN") == -1: # some reaction are classified as pathways
          if i == 0:
            pathways.append(str(pw.common_name))
          else:
            pathways.append(str(pw.common_name))
          i += 1
    else:
      if str(pws).find("RXN") == -1:
        pathways.append(no_style(str(pws.common_name)))
  return pathways


def metabolite_name(org, metabolite):
  """Returns full name of a metabolite"""
  return org.get_name_string(metabolite)


def metabolite_formula(org, metabolite):
  if "CHEMICAL-FORMULA" in metabolite.keys():
    formula = str(metabolite.chemical_formula)
    if formula == "None": 
      formula=""
    else: 
      formula = formula.replace("[","").replace("]","").replace("'","").replace(", ","")
    return formula
  else: 
    return ""


def metabolite_compartment(org, reaction, metabolite, side):
  """Returns metabolite's compartment given a reaction"""
  #sm.transporter(reaction) is not working allways!!
  compartments = org.compartments_of_reaction(reaction)
  if len(compartments) == 1:  # <=> is not a transporter
    return compartments[0]
  else:                       # <=> is a transporter
    val = org.get_value_annot(reaction, side, metabolite, "compartment") 
    if val == None: # sadly happens sometimes...
      print "database has no compartment info for", metabolite, "in reaction", reaction
      comp = "CCO-IN" # assuming cytosol
      print "... assuming", comp
      return comp
    else:
      return str(val)


def reaction_meta_stoich(org, reaction):
  """Anlayses a given reaction and returns its metabolites with stoichiometry"""
  meta_stoich_dic = {}
  reactants = org.reaction_reactants_and_products(reaction)[0]
  for reactant in reactants:
    stoich = org.get_value_annot(reaction, "left", reactant, "coefficient") # stoichiometry
    if stoich == None: stoich = 1 # None <=> 1
    compartment = metabolite_compartment(org, reaction, reactant, "left")
    formula = metabolite_formula(org, reactant)
    name = metabolite_name(org, reactant)
    metabolite = Metabolite(str(reactant), formula, name, compartment)
    meta_stoich_dic[metabolite] = -int(stoich) # negative because reactant is consumed
  products = org.reaction_reactants_and_products(reaction)[1]
  for product in products:
    stoich = org.get_value_annot(reaction, "right", product, "coefficient") # stoichiometry
    if stoich == None: stoich = 1 # None <=> 1
    compartment = metabolite_compartment(org, reaction, product, "right")
    formula = metabolite_formula(org, product)
    name = metabolite_name(org, product)
    metabolite = Metabolite(str(product), formula, name, compartment)
    meta_stoich_dic[metabolite] = int(stoich)
  return meta_stoich_dic


def reaction_reversible(org, reaction):
  """Returns True/False if reaction is reversible or not"""
  dir = reaction.reaction_direction
  if dir == None: # metacyc has no entry if it's reversible
    return True
  else:
    return False


def reaction_gene_reaction_rule(org, reaction):
  """Returns logical link between reaction's involved genes
     and <=> enzyme comples, or <=> isoenzyme"""
  gpr = ""
  i = 0
  re_enzymes = org.enzymes_of_reaction(reaction)
  if re_enzymes != None:
    if hasattr(re_enzymes, '__iter__'): # check if iterable
      for enzyme in re_enzymes:
        if org.complex(enzyme) != True: # if it's not a protein complex -> isoenzymes (or)
          if i == 0:
            gpr   += str(org.genes_of_protein(enzyme))
          else:
            gpr   += " or " + str(org.genes_of_protein(enzyme))
        else: # if it's a protein complex -> no isoenzymes (and)
          complex = org.genes_of_protein(enzyme)
          if hasattr(complex, '__iter__'): # check if iterable
            for subunit in complex:
              if subunit == complex[0]: # first in list -> brackets
                if i == 0:
                  gpr   += "(" + str(subunit)
                else:
                  gpr   += " or (" + str(subunit)
              else:
                gpr   += " and " + str(subunit)
              if subunit == complex[len(complex) - 1]: #last in list -> brackets
                gpr   += ")"
          else: # complex with no subunit (gene not avaible)
            if i == 0:
              gpr   += "NA"
            else:
              gpr   += " or " + "NA"
        i+=1
    else:
      gpr   = str(org.genes_of_protein(re_enzymes)[0])
      
  gpr   = re.sub("[\[\]]", "", gpr) # remove [ ] brackets
  return "("+gpr+")"

