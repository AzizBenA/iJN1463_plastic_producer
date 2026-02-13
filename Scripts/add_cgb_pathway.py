from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction, Metabolite
from pathlib import Path

if __name__ == '__main__':
    polymer_length = 1
    model = Path('Models/250709_iJN1463_updated_5.sbml')
    cph_production_rxn = Reaction(
        'CPHA1', name = 'Cyanophycin Synthetase from Anabaena sp. PCC 7120'
    )
    cph_production_rxn.gene_reaction_rule = 'all3879'
    cyanophycin = Metabolite('cgp_c',name= 'cyanophycin', compartment='c')
    cph_production_rxn.add_metabolites(
        {model.metabolites.asp__L_c: -(polymer_length+1),
         model.metabolites.arg__L_c: -polymer_length,
         model.metabolites.ATP_c: -polymer_length,
         model.metabolites.ADP_c: polymer_length,
         model.metabolites.pi_c: polymer_length,
         model.metabolites.h_c: polymer_length,
         cyanophycin: 1})
    model.add_reactions([cph_production_rxn])
    model.add_boundary(cyanophycin, type = 'sink')