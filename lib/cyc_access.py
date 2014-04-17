import re
from cobra import Model, Reaction, Metabolite
import itertools
import copy

compartment_dic = {"CCO-IN":"c", "CCO-EXTRACELLULAR":"e", "CCO-CYTOSOL":"c", "CCO-PERI-BAC":"p", "CCO-OUT":"e"}



def test():
  print "hello"


def reaction_name(org, reaction):
  """Returns the full name of a reaction"""
  name = reaction.common_name
  if name == None: name = reaction.systematic_name
  if name == None:
    enzyme = reaction.enzymatic_reaction
    if isinstance(enzyme, list):
      name = enzyme[0].common_name 
    elif enzyme != None:
      name = enzyme.common_name 
  if name == None: name = str(reaction)
  return no_style(name)


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


def id_cleaner(string):
  """Returns a string without -.+ (special charakters in ptools)"""
  return string.replace(".","").replace("-","").replace("+","").replace("|","")


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
  return no_style(org.get_name_string(metabolite))


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
    return compartment_dic[str(compartments[0])]
  else:                       # <=> is a transporter
    val = org.get_value_annot(reaction, side, metabolite, "compartment") 
    if val == None: # sadly happens sometimes...
      print "database has no compartment info for", metabolite, "in reaction", reaction
      comp = "CCO-IN" # assuming cytosol
      print "... assuming", comp
      return compartment_dic[comp]
    else:
      return compartment_dic[str(val)]


def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False


def reaction_meta_stoich(org, reaction, substitutions):
  """Anlayses a given reaction and returns its metabolites with stoichiometry"""
  meta_stoich_dic = {}
  if reaction.reaction_direction == "RIGHT-TO-LEFT": # Attention: pathwaytools considers reactions  of the type "B <- A" confusingly (B is reactant and A is product!!)
    reactants = org.reaction_reactants_and_products(reaction)[1]
    products = org.reaction_reactants_and_products(reaction)[0]
  else:
    reactants = org.reaction_reactants_and_products(reaction)[0]
    products = org.reaction_reactants_and_products(reaction)[1]

  for reactant in reactants:
    if id_cleaner(str(reactant)) in substitutions.keys(): reactant = org.get_frame_labeled(substitutions[id_cleaner(str(reactant))])[0] # get substitution of a metabolite if avaible
    stoich = org.get_value_annot(reaction, "left", reactant, "coefficient") # stoichiometry
    if stoich == None: stoich = 1 # None <=> 1
    compartment = metabolite_compartment(org, reaction, reactant, "left")
    formula = metabolite_formula(org, reactant)
    name = metabolite_name(org, reactant)
    abbr = id_cleaner(str(reactant)+"_"+compartment)
    metabolite = Metabolite(abbr, formula, name, compartment)
    if is_number(stoich):
      meta_stoich_dic[metabolite] = -int(stoich) # negative because reactant is consumed
    else:
      meta_stoich_dic[metabolite] = stoich
  for product in products:
    if id_cleaner(str(product)) in substitutions.keys(): product = org.get_frame_labeled(substitutions[id_cleaner(str(product))])[0] # get substitution of a metabolite if avaible
    stoich = org.get_value_annot(reaction, "right", product, "coefficient") # stoichiometry
    if stoich == None: stoich = 1 # None <=> 1
    compartment = metabolite_compartment(org, reaction, product, "right")
    formula = metabolite_formula(org, product)
    name = metabolite_name(org, product)
    abbr = id_cleaner(str(product)+"_"+compartment)
    metabolite = Metabolite(abbr, formula, name, compartment)
    if is_number(stoich): 
      meta_stoich_dic[metabolite] = int(stoich)
    else:
      meta_stoich_dic[metabolite] = stoich
  return meta_stoich_dic


def reaction_reversible(org, reaction):
  """Returns True/False if reaction is reversible or not"""
  dir = reaction.reaction_direction
  if str(dir).lower() == "reversible" or dir == None: # consider reaction with unknow reversibility as reversible
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
  #return "("+gpr+")"
  return ("( "+gpr+" )").replace("-","")


def reaction_is_generic(org, reaction, exceptions, substitutions):
  """Returns True if a reaction contains at least one generic/unspecific metabolite"""
  all_metabolites = org.reaction_reactants_and_products(reaction)[0] + org.reaction_reactants_and_products(reaction)[1]
  for metabolite in all_metabolites:
    if id_cleaner(str(metabolite)) in substitutions.keys(): 
      metabolite = org.get_frame_labeled(substitutions[id_cleaner(str(metabolite))])[0] # get substitution of a metabolite if avaible
    if id_cleaner(str(metabolite)) not in exceptions and org.is_class(metabolite):
      return True
  return False


def reaction_get_generic(org, reaction, exceptions, substitutions):
  """Returns a set with generic/unspecific metabolites of a reaction"""
  generic = set()
  all_metabolites = org.reaction_reactants_and_products(reaction)[0] + org.reaction_reactants_and_products(reaction)[1]
  for metabolite in all_metabolites:
    if id_cleaner(str(metabolite)) in substitutions.keys(): 
      metabolite = org.get_frame_labeled(substitutions[id_cleaner(str(metabolite))])[0] # get substitution of a metabolite if avaible
    if id_cleaner(str(metabolite)) not in exceptions and org.is_class(metabolite):
      generic.add(str(metabolite))
  return generic


def find_specific(org, generic_metabolite):
  """Returns a list of specific metabolites given a generic one
  generic metabolites are classes (eg. fatty acids) which should be avoided because they are ambiguous"""
  specified = []
  if org.is_class(generic_metabolite): # if it's a class
    if hasattr(org.get_class_all_instances(generic_metabolite), '__iter__'): # if this class has instances
      for c in org.get_class_all_instances(generic_metabolite):
        if c not in specified: specified.append(c)
    if org.get_class_all_subs(generic_metabolite) != None: # if this class has subclases
      for sub in org.get_class_all_subs(generic_metabolite):
        for f in find_specific(org, sub):
          if f not in specified: specified.append(f)
  return specified


def meta_stoich_replace(dic, old, new):
  """Returns a dictionary in which old key is exchanged with new key"""
  dic_new = {}
  for member in dic:
    if member == old: dic_new[new] = dic[old]
    else: dic_new[member] = dic[member]
  return dic_new


def metabolite_from_string(setlistdic, string):
  """Returns an object from class Metabolite if found in a set, list, dic, otherwise None"""
  for entry in setlistdic:
    if entry.id == string:
      return entry
    if entry.id == string + "_" + entry.compartment: # take care of ids which have a compartment tag
      return entry
  return None


def reaction_generic_specified(org, reaction, org_reaction, generic_exceptions, substitutions, generic_assignment):
  specified_reactions     = []  # list of new specified reactions
  all_metabolites         = org.reaction_reactants_and_products(reaction)[0] + org.reaction_reactants_and_products(reaction)[1]
  generics_substitutions  = {}  # dictionary which contains a list with specific metabolites for each generic one
  list_generics           = []  # list of generic metabolite 
  multilist_specifics     = []  # list containing a list for every generic metabolite containing its specific metabolites
  meta_stoich             = reaction_meta_stoich(org, reaction, substitutions)
  for metabolite in all_metabolites:
    if str(metabolite).replace("|","") in substitutions.keys(): metabolite = org.get_frame_labeled(substitutions[str(metabolite).replace("|","")])[0] # get substitution of a metabolite if avaible
    if org.is_class(metabolite) and id_cleaner(str(metabolite)) not in list_generics and id_cleaner(str(metabolite)) not in generic_exceptions:
      specifics = find_specific(org, metabolite)
      generics_substitutions[id_cleaner(str(metabolite))] = specifics # cleaning troubling characters in metacyc id names
      list_generics.append(id_cleaner(str(metabolite)))
      multilist_specifics.append(specifics)
 
  tmp = {id_cleaner(value): len(generics_substitutions[value]) for value in generics_substitutions if len(generics_substitutions[value]) > 1} # only abstract metabolites with more than 1 specifification
  tmp_set = set() # specified set of the complex case
  for value in generics_substitutions:
    if len(generics_substitutions[value]) > 1:
      tmp_set |= set(map(str, generics_substitutions[value]))
  #print len(generics_substitutions)
  
  combinations = itertools.product(*multilist_specifics)  # all combinations of specifications in a reaction (k-combination, no order, without replacement)
  
  if len(tmp) > 1: # complex reaction possbible due to some equivocal assignments (e.g. NAD(P)-NAD(P)H => specification procedure has to map NAD-NADH and NADP-NADPH, not NAD-NADPH and not NADP-NADH)
    # try to solve complex reaction with equivocal assignment by reading file with common mappings for each case
    print "generics:", tmp
    print "identified specifications:",tmp_set
    if tmp_set <= set(generic_assignment.keys()): # check if for all complex parts of reaction there are generic assignments => specification is possible
      combinations_new = []
      for c in combinations: # for each combination 
        #print c
        found = False
        for element in map(str, c): # for each element of the combination
          if element in tmp_set: # if element is in specified set of the complex case
            if not generic_assignment[element] in map(str, c):# check if assignment is complete (i.e. if A is present the corresponding B has to be present
              found = False
              break
            else: 
              found = True
        if found: 
          combinations_new.append(c)
      print "useful pairs:", combinations_new
      combinations = combinations_new
    else:
      print "Complex reaction:", org_reaction, org_reaction.reaction, "\nnot added!"
      print tmp
      return []
  nr = 0
  for combination in combinations:
    nr += 1
    new_meta_stoich = {}
    for index, generic in enumerate(list_generics):
      specific = combination[index]
      for entry in meta_stoich: # change generic to specific metabolite in reaction list (meta_stoich) and build a new reaction
        if entry.id[:entry.id.find("_")] == str(generic).replace("|",""):
          generic_metabolite      = entry
          specific_metabolite     = copy.deepcopy(generic_metabolite)
          #specific_metabolite.id  = specific_metabolite.id.replace(str(generic).replace("|",""), str(specific))
          specific_metabolite.id  = specific_metabolite.id.replace(generic, str(specific))
          specific_metabolite.name= metabolite_name(org, specific)
          new_meta_stoich[specific_metabolite] = meta_stoich[entry]
        elif entry.id[:entry.id.find("_")] not in list_generics:
          new_meta_stoich[entry] = meta_stoich[entry]
    reaction_new = Reaction(org_reaction.id + "_" + str(nr))
    reaction_new.name                   = org_reaction.name
    reaction_new.subsystem              = org_reaction.subsystem 
    reaction_new.lower_bound            = org_reaction.lower_bound
    reaction_new.upper_bound            = org_reaction.upper_bound
    reaction_new.objective_coefficient  = org_reaction.objective_coefficient
    reaction_new.add_metabolites(new_meta_stoich)
    reaction_new.gene_reaction_rule     = org_reaction.gene_reaction_rule
    #print reaction_new.reaction
    specified_reactions.append(reaction_new)
  return specified_reactions


def read_generic_exceptions(filename):
  """returns a list with possible generic metabolites which shouldnt be specified"""
  list = []
  file = open(filename, "r")
  for line in file:
    if line != "\n" and line.lstrip()[0] != "#":
      name = line.rstrip("\n")
      print "added exception for generic metabolite", name
      list.append(id_cleaner(name))
  return list


def substitutions_dic(filename):
  """returns a dictionary containing old:new substitutions of metabolites"""
  dic = {}
  file = open(filename, "r")
  for line in file:
    if line != "\n" and line.lstrip()[0] != "#":
      split = line.rstrip("\n").split(":")
      if len(split) == 2:
        old = split[0]
        new = split[1]
        print old, "is going to be substituted with", new
        dic[id_cleaner(old)]=new # attention: only the key of the dictionary can be cleaned (no annoying chars) because the value is used to access ptools objects via api
      else: print filename, "error in line", line
  return dic


def get_diffusion_reactions(org, filename):
  """returns a list of diffusion reaction read from file"""
  reaction_list = []
  file = open(filename, "r")
  for line in file:
    if line != "\n" and line.lstrip()[0] != "#":
      split = line.rstrip("\n").split(":")
      if len(split) == 2 and compartment_dic.has_key(split[0]):
        compartment = compartment_dic[split[0]]
        substance  = split[1]
        print "added diffusion reaction for", substance, "into", compartment 
        metacyc_metabolite = org.get_frame_labeled(substance)[0]
        formula = metabolite_formula(org, metacyc_metabolite)
        name = metabolite_name(org, metacyc_metabolite)
        abbr = id_cleaner(str(metacyc_metabolite)+"_"+compartment)
        metabolite = Metabolite(abbr, formula, name, compartment)
        reaction = Reaction("EX_"+str(metabolite))
        reaction.name = "Diffusion of " + name
        reaction.lower_bound  = -1000 
        reaction.upper_bound  = 1000
        meta_stoich = {metabolite:-1}
        reaction.add_metabolites(meta_stoich)
        reaction_list.append(reaction)
      else: print filename, "error in line", line
  return reaction_list


def get_generic_assignment(filename):
  """returns a dictionary containing for each generic compounds (for which assignments exists) a tuple with its assignments"""
  dic_assignments = {}
  file = open(filename, "r")
  for line in file:
    if line != "\n" and line.lstrip()[0] != "#":
      split = line.rstrip("\n").split(":")
      #if len(split) == 3:
      if len(split) == 2:
        #generic_name  = id_cleaner(split[0])  # name of generic compound
        #assignment1   = id_cleaner(split[1])  # first metabolite to be assigned
        assignment1   = id_cleaner(split[0])  # first metabolite to be assigned
        #assignment2   = id_cleaner(split[2])  # second metabolite to be assigned
        assignment2   = id_cleaner(split[1])  # second metabolite to be assigned
        #dic_assignments[generic_name] = assignment1, assignment2
        dic_assignments[assignment1] = assignment2
        dic_assignments[assignment2] = assignment1
      else: print filename, "error in line", line
  return dic_assignments
