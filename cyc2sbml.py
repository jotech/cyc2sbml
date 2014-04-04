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
import operator

r_ignored = open("reactions_ignored", "w")
p_ignored = open("pathways_ignored", "w")
m_generic = open("harmul_generics", "w")
p_ignored_set = set() # set of ignored pathways
r_ignored_set = set() # set of ignored reactions
m_generic_set = {} # dictionary of generic metabolites and how often they lead to a irgnored reaction

print pycyc.all_orgs()

answer_org  = raw_input("Which Organism shoulb be exportet to sbml? ")
org         = pycyc.open(answer_org)
print "Loading", answer_org, "with", len(org.all_rxns(":all")), "reactions"

answer_generic = raw_input("\nDo you want to take care of generic metabolites (e.g. tRNA, carboxylates)?\n All subclasses of a generic metabolite will be searched for instances (that is specific metabolites) and specific reactions (could be much more!) will be added instead of the generic one\n [y/n] ")
care_generics = True if answer_generic == "y" else False
if care_generics:
  answer_exceptions = raw_input("\nAre some metabolites (./exceptions.txt) to be kept even if they are generic metabolites? This could be usefull e.g. to summarize lipid metabolism.\n [y/n] ")
  generic_exceptions = cyc.read_generic_exceptions("./exceptions.txt") if answer_exceptions == "y" else []
else: generic_exceptions = []
answer_substitutions = raw_input("\nSubstitutions defined in ./substitutions.txt could be read and applied to exchange certain metabolites in pathwaytools database. Is this to be done? [y/n] ")
substitutions = cyc.substitutions_dic("./substitutions.txt") if answer_substitutions == "y" else {}
answer_diffusion = raw_input("\nShould a exchange reaction for membrane permeable substances (defined in ./diffusion.txt) be added automatically to the model? [y/n] ")
if answer_diffusion == "y":
  diffusion_reactions_list = cyc.get_diffusion_reactions(org, "./diffusion.txt")
  diffusion_reactions = True
else: diffusion_reaction = False

answer_start = raw_input("\n---\nReady to start? [y/n] ")
if not answer_start == "y": quit()

model = Model(answer_org)

#for r in org.all_rxns(":all"): # all reaction
#for r in [org.get_frame_labeled("R601-RXN")[0]]:  # consider only some reaction for testing
for r in org.all_rxns(":metab-smm") + org.all_rxns(":transport"): # only metabolic reactions 
#for r in org.all_rxns(":all")[0:10]: # only the first reactions -> debugging
  reaction                        = Reaction(cyc.id_cleaner(str(r)))
  reaction.name                   = cyc.reaction_name(org, r)
  reaction.subsystem              = cyc.reaction_subsystem(org, r)
  reaction.lower_bound            = -1000 if cyc.reaction_reversible(org, r) else 0
  reaction.upper_bound            = 1000
  reaction.objective_coefficient  = 0
  reaction.add_metabolites(cyc.reaction_meta_stoich(org, r, substitutions))
  reaction.gene_reaction_rule     = cyc.reaction_gene_reaction_rule(org, r)

  #print reaction.print_values()
  #if reaction.check_mass_balance() != []: print reaction, "is not balanced"  # throws error sometimes!

  if care_generics and cyc.reaction_is_generic(org, r, generic_exceptions, substitutions):
    specific_reactions = cyc.reaction_generic_specified(org, r, reaction, generic_exceptions, substitutions)
    model.add_reactions(specific_reactions)
    print "\nabstract reaction:", reaction, reaction.reaction, "\n\tadded", len(specific_reactions), "specific reactions"
    if len(specific_reactions) == 0:
      r_ignored_set.add(str(r))
      print >>r_ignored, str(r), reaction.reaction
      if r.in_pathway != None:
        plist = r.in_pathway if isinstance(r.in_pathway, list) else [r.in_pathway]
        p_ignored_set |= set((map(str, plist)))
      for g in cyc.reaction_get_generic(org, r, generic_exceptions, substitutions):
        m_generic_set[g] = m_generic_set.get(g, 0) + 1 # count for each generic metabolite how often it causes a irgnored reaction
  else:
    model.add_reaction(reaction)

if diffusion_reactions: model.add_reactions(diffusion_reactions_list) # adding automatically additional diffusion reactions

for pwy in p_ignored_set: print >>p_ignored, pwy, org.get_name_string(pwy) 

for s in sorted(m_generic_set.iteritems(), key=operator.itemgetter(1)): print >> m_generic, s[0], s[1]

print '\n---\n%i reaction in model' % len(model.reactions)
print '%i metabolites in model' % len(model.metabolites)
print '%i genes in model\n---\n' % len(model.genes)

print "ignored reactions:", len(r_ignored_set)
print "thus incomplete pathways:", len(p_ignored_set), "\n"

sbml_out_file = answer_org+'.xml'
sbml_level    = 2
sbml_version  = 1  # Writing level 2, version 4 is not completely supported.
write_cobra_model_to_sbml_file(model, sbml_out_file, sbml_level, sbml_version, print_time=True)
