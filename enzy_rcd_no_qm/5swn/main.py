import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh

prefix='../../starting_points/5swn/'

sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/full_chain_renumbered.pdb')
#
ADP = lp.get_ligand(f"{prefix}/FAH.mol2")
DAR = lp.get_ligand(f"{prefix}/DAR_protonated.mol2")
NO3 = lp.get_ligand(f"{prefix}/NO3.mol2")
MG  = eh.get_metal('MG', charge=2)
#
stru.add(FAH, chain_name='Z', net_charge= -1, multiplicity=1) 

help(eh.preparation.place_ligand)
exit( 0 )
eh.preparation.place_ligand(stru, 'Z', 1)
#eh.preparation.place_ligand(stru, 'B', 1)
#
#sp.save_structure('temp.pdb', stru)
#stru = sp.get_structure('temp.pdb')

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")

eh.dock_reactants( stru,
                    [stru.get('Y.1'), stru.get('Z.1')],
                    constraints=constraints,
                    n_struct=100,
                    cluster_distance=2.5,
                    clash_cutoff=3,
                    cst_energy=2000,
                    )
