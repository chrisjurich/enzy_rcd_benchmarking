import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh


#session = eh.interface.pymol.new_session()
#
#eh.interface.pymol.general_cmd(session, [('fetch', 'FAH'), ('save', 'FAH.mol2')])
#
#lp = eh.Mol2Parser()
#
#eh.interface.moe.protonate('FAH.mol2')
#exit( 0 )
#prefix='../../starting_points/1bg0/'
#
sp = eh.PDBParser()

stru = sp.get_structure('3R3V.pdb')
print(stru)
eh.preparation.add_missing_residues(stru, eh.preparation.identify_missing_residues('3R3V'))
print(stru)
exit( 0 )
h.preparation.protonate_stru(stru)

sp.save_structure('full_chain_renumbered.pdb', stru)




print(stru)
#lp = eh.Mol2Parser()
#
#stru = sp.get_structure(f'{prefix}/full_chain_renumbered.pdb')
##
#ADP = lp.get_ligand(f"{prefix}/ADP.mol2")
#DAR = lp.get_ligand(f"{prefix}/DAR_protonated.mol2")
#NO3 = lp.get_ligand(f"{prefix}/NO3.mol2")
#MG  = eh.get_metal('MG', charge=2)
##
#stru.add(MG , chain_name="B", net_charge= 2, multiplicity=1)
#stru.add(ADP, chain_name='X', net_charge=-3, multiplicity=1)
#stru.add(NO3, chain_name='Y', net_charge=-1, multiplicity=1) 
#stru.add(DAR, chain_name='Z', net_charge= 1, multiplicity=1) 
#
##eh.preparation.place_ligand(stru, 'X', 1)
##eh.preparation.place_ligand(stru, 'B', 1)
##
##sp.save_structure('temp.pdb', stru)
##stru = sp.get_structure('temp.pdb')
#
#constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")
#
#eh.dock_reactants( stru,
#                    [stru.get('Y.1'), stru.get('Z.1')],
#                    constraints=constraints,
#                    n_struct=100,
#                    cluster_distance=2.5,
#                    clash_cutoff=3,
#                    cst_energy=2000,
#                    )
