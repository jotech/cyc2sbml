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
  reaction                        = Reaction(cyc.id_cleaner(str(r)))
  reaction.name                   = cyc.reaction_name(org, r)
  reaction.subsystem              = cyc.reaction_subsystem(org, r)
  reaction.lower_bound            = -1000 if cyc.reaction_reversible(org, r) else 0
  reaction.upper_bound            = 1000
  reaction.objective_coefficient  = 0
  reaction.add_metabolites(cyc.reaction_meta_stoich(org, r))
  reaction.gene_reaction_rule     = cyc.reaction_gene_reaction_rule(org, r)

  #print reaction.print_values()
  #if reaction.check_mass_balance() != []: print reaction, "is not balanced"  # throws error sometimes!

  if care_generics and cyc.reaction_is_generic(org, r):
    specific_reactions = cyc.reaction_generic_specified(org, r, reaction)
    model.add_reactions(specific_reactions)
    print "\nabstract reaction:", reaction, reaction.reaction, "\n\tadded", len(specific_reactions), "specific reactions"
    list = cyc.reaction_generic_specified(org, r, reaction)
    for l in list:
      print l, l.reaction
    #import pdb; pdb.set_trace()
  else:
    model.add_reaction(reaction)

print '\n---\n%i reaction in model' % len(model.reactions)
print '%i metabolites in model' % len(model.metabolites)
print '%i genes in model\n---\n' % len(model.genes)

sbml_out_file = answer_org+'.xml'
sbml_level    = 2
sbml_version  = 1  # Writing level 2, version 4 is not completely supported.
write_cobra_model_to_sbml_file(model, sbml_out_file, sbml_level, sbml_version, print_time=True)
