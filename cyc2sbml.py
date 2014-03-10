# cyc2sbml.py
#
# This program reads out an organism's database from pathwaytools (e.g. metacyc) and writes a sbml model file
#
# - Assuming that ptools is running in api mode: ./pathway-tools -lisp -api
# - Needs external python packages: pycyc, cobrapy


import pycyc
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import write_cobra_model_to_sbml_file
import cyc_access as cyc


print pycyc.all_orgs()

answer_org  = raw_input("Which Organism shoulb be exportet to sbml? ")
org         = pycyc.open(answer_org)
print "Loading", answer_org, "with", len(org.all_rxns(":all")), "reactions"

answer_generic = raw_input("\nDo you want to take care of generic metabolites (e.g. tRNA, carboxylates)?\n All subclasses of a generic metabolite will be searched for instances (that is specific metabolites) and specific reactions (could be much more!) will be added instead of the generic one\n [y/n] ")
care_generics = True if answer_generic == "y" else False

model = Model(answer_org)

#for r in org.all_rxns(":all"):
for r in org.all_rxns(":all")[0:10]: # only the first reactions -> debugging
  reaction                        = Reaction(str(r))
  reaction.name                   = cyc.reaction_name(org, r)
  reaction.subsystem              = cyc.reaction_subsystem(org, r)
  reaction.lower_bound            = 0
  reaction.upper_bound            = 1000
  reaction.objective_coefficient  = 0
  reaction.reversibility          = cyc.reaction_reversible(org, r)
  reaction.add_metabolites(cyc.reaction_meta_stoich(org, r))
  reaction.add_gene_reaction_rule(cyc.reaction_gene_reaction_rule(org, r))

  if care_generics and cyc.reaction_is_generic(org, r):
    print "abstract reaction:", reaction.reaction
    #model.add_reactions(cyc.reaction_generic_specified(org, r))
    list = cyc.reaction_generic_specified(org, r, reaction)
    for l in list:
      print l.reaction
  else:
    model.add_reactions(reaction)


print '%i reaction in model' % len(model.reactions)
print '%i metabolites in model' % len(model.metabolites)
print '%i genes in model' % len(model.genes)

sbml_out_file = answer_org+'.xml'
sbml_level    = 3
sbml_version  = 1  # Writing level 1, version 4 is not completely supported.
#write_cobra_model_to_sbml_file(model, sbml_out_file, sbml_level, sbml_version, print_time=True)
