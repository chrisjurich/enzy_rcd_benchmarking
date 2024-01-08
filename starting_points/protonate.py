import enzy_htp as eh




parser = eh.PDBParser()


stru = parser.get_structure('./full_chain_renumbered.pdb')

eh.preparation.protonate_stru(stru)

parser.save_structure('start.pdb', stru)

