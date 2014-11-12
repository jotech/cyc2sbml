# cyc2sbml.py
#
# This program reads out an organism's database from pathwaytools (e.g. metacyc) and writes a sbml model file
#
# - Assuming that ptools is running in api mode: ./pathway-tools -lisp -api
# - Needs external python packages: pycyc, cobrapy


import pycyc
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model
import sys
import os
sys.path.append(os.path.abspath("./lib"))
import cyc_access as cyc
import operator

r_ignored = open("reactions_ignored", "w")
p_ignored = open("pathways_ignored", "w")
m_generic = open("harmul_generics", "w")
mass_balance = open("mass_balance", "w")
p_ignored_set = {} # dictionary of ignored pathways and blocking metabolites
r_ignored_set = set() # set of ignored reactions
m_generic_set = {} # dictionary of generic metabolites and how often they lead to a ignored reaction
p_generic_set = {} # dictionary of ignored pathways and how much ignored generic metabolites are responsible
r_generic = 0 # count generic reactions
r_total   = 0 # count reactions


#
# I. input
#

print pycyc.all_orgs()

answer_org  = raw_input("Which Organism shoulb be exportet to sbml? ")
org         = pycyc.open(answer_org)
print "Loading", answer_org, "with", len(org.all_rxns(":all")), "reactions"

answer_generic = raw_input("\nDo you want to take care of generic metabolites (e.g. tRNA, carboxylates)?\n All subclasses of a generic metabolite will be searched for instances (that is specific metabolites) and specific reactions (could be much more!) will be added instead of the generic one\n [y/n] ")
care_generics = True if answer_generic == "y" else False
if care_generics:
  answer_exceptions = raw_input("\nAre some metabolites (./conf/exceptions.txt) to be kept even if they are generic metabolites? This could be usefull e.g. to summarize lipid metabolism.\n [y/n] ")
  generic_exceptions = cyc.read_generic_exceptions("./conf/exceptions.txt") if answer_exceptions == "y" else []
else: generic_exceptions = []

answer_substitutions = raw_input("\nSubstitutions defined in ./conf/substitutions.txt could be read and applied to exchange certain metabolites in pathwaytools database. Is this to be done? [y/n] ")
substitutions = cyc.substitutions_dic("./conf/substitutions.txt") if answer_substitutions == "y" else {}

answer_diffusion = raw_input("\nShould a exchange reaction for membrane permeable substances (defined in ./conf/diffusion.txt) be added automatically to the model? [y/n] ")
if answer_diffusion == "y":
  diffusion_reactions_list = cyc.get_diffusion_reactions(org, "./conf/diffusion.txt")
  diffusion_reactions = True
else: diffusion_reactions = False

answer_bigg_names = raw_input("\nDo you want to use bigg reaction names? [y/n] ")
bigg_names = True if answer_bigg_names == "y" else False
if bigg_names: bigg_reaction_dic = cyc.get_bigg_reaction_dic("./conf/metacyc_bigg.txt")

answer_bigg_names_metabolites = raw_input("\nDo you want to use bigg metabolite names? [y/n] ")
bigg_names_metabolites = True if answer_bigg_names_metabolites == "y" else False
metabolites_dic = cyc.get_bigg_metabolites_dic("./conf/metacyc_bigg_substances.txt") if bigg_names_metabolites else {}

answer_ignore = raw_input("\nShould some reactions defined in ./conf/ignore.txt be ignored? [y/n] ")
to_ignore = cyc.get_to_ignore_reactions("./conf/ignore.txt") if answer_ignore == "y" else set()

answer_gene_names = raw_input("\nDo you want to use different gene names (locus tags) as used in ptools? .conf/gene_names.txt [y/n] ")
gene_names = cyc.get_gene_names_dic("./conf/gene_names_dict.txt") if answer_gene_names == "y" else {}

#
# II. making of
#

answer_start = raw_input("\n---\nReady to start? [y/n] ")
if not answer_start == "y": quit()
else: print "\n\n"

model = Model(answer_org)

#for r in org.all_rxns(":all"): # all reaction
#for r in [org.get_frame_labeled("DNA-DIRECTED-DNA-POLYMERASE-RXN")[0]]:  # consider only some reaction for testing
for r in org.all_rxns(":metab-all") + org.all_rxns(":transport"): # only metabolic reactions 
#for r in org.all_rxns(":all")[0:10]: # only the first reactions -> debugging
  if str(r) in to_ignore: continue # if reaction is defined to be ignored 
  r_total += 1
  if bigg_names and bigg_reaction_dic.has_key(str(r)):
    reaction                      = Reaction(bigg_reaction_dic[str(r)]) if bigg_reaction_dic[str(r)] not in model.reactions else Reaction(cyc.id_cleaner(str(r)))
  else:
    reaction                      = Reaction(cyc.id_cleaner(str(r)))
  if reaction.id == "": reaction.id = r.ec_number
  reaction.name                   = cyc.reaction_name(org, r) + " " + str(r.ec_number)
  reaction.subsystem              = cyc.reaction_subsystem(org, r)
  reaction.lower_bound            = -1000 if cyc.reaction_reversible(org, r) else 0
  reaction.upper_bound            = 1000
  reaction.objective_coefficient  = 0
  reaction.add_metabolites(cyc.reaction_meta_stoich(org, r, substitutions))
  reaction.gene_reaction_rule     = cyc.reaction_gene_reaction_rule(org, r, gene_names)

  #print reaction.print_values()
  #if reaction.check_mass_balance() != []: print reaction, "is not balanced"  # throws error sometimes!

  if care_generics and cyc.reaction_is_generic(org, r, generic_exceptions, substitutions):
    r_generic += 1
    specific_reactions = cyc.reaction_generic_specified(org, r, reaction, generic_exceptions, substitutions, cyc.get_generic_assignment("./conf/generic_assignment.txt"))
    model.add_reactions(specific_reactions)
    #print "\nabstract reaction:", str(r), reaction.reaction, "\n\tadded", len(specific_reactions), "specific reactions"
    print "\nabstract reaction:", str(r),"\n\tadded", len(specific_reactions), "specific reactions"
    if len(specific_reactions) == 0:
      r_ignored_set.add(str(r))
      #print >>r_ignored, str(r), reaction.reaction
      generics = cyc.reaction_get_generic(org, r, generic_exceptions, substitutions)
      if r.in_pathway != None: # remember missed pathways
        plist = r.in_pathway if isinstance(r.in_pathway, list) else [r.in_pathway]
        for path in plist:
          #p_ignored_set |= set((map(str, plist)))
          if p_ignored_set.has_key(str(path)): p_ignored_set[str(path)] = p_ignored_set[str(path)] | generics
          else: p_ignored_set[str(path)] = generics
      for g in generics: # remember missed reactions
        m_generic_set[g] = m_generic_set.get(g, 0) + 1 # count for each generic metabolite how often it causes a irgnored reaction
  else:
    model.add_reaction(reaction)

if diffusion_reactions: model.add_reactions(diffusion_reactions_list) # adding automatically additional diffusion reactions

if bigg_names_metabolites: model = cyc.change_metabolite_names(model, metabolites_dic)

#
# III. output
# 

for r in model.reactions: 
  if r.check_mass_balance() != []: 
    print >>mass_balance, r.id, r.name, r.check_mass_balance()

for pwy in p_ignored_set: 
  print >>p_ignored, pwy, org.get_name_string(pwy), p_ignored_set[pwy]
  p_generic_set[org.get_name_string(pwy)] = len(p_ignored_set[pwy])

for s in sorted(m_generic_set.iteritems(), key=operator.itemgetter(1)): print >> m_generic, s[0], s[1]
for s in sorted(p_generic_set.iteritems(), key=operator.itemgetter(1)): print >> m_generic, s[0], s[1]

print '\n---\n%i reaction in model' % len(model.reactions)
print '%i metabolites in model' % len(model.metabolites)
print '%i genes in model\n---\n' % len(model.genes)

print "reactions in database:", r_total
print "generic reactions:", r_generic
print "generic metabolites:", len(m_generic_set.keys())
print "ignored reactions:", len(r_ignored_set)
print "thus incomplete pathways:", len(p_ignored_set), "\n"

sbml_out_file = answer_org+'.xml'
write_sbml_model(model, sbml_out_file, use_fbc_package=False)
