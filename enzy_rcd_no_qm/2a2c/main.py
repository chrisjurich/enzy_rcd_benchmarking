import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh


eh.config['mole2.MONO'] = "/accre/arch/easybuild/software/Compiler/GCCcore/6.4.0/Mono/5.4.1.6/bin/mono"
eh.config['mole2.PROBE'] = 5.0
eh.config['mole2.INNER'] = 1.05
eh.config['mole2.IGNORE_HETATM'] = True


prefix='../../starting_points/2a2c/'



sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/full_chain_renumbered.pdb')
ADP = lp.get_ligand(f"./ADP.mol2")
NG1 = lp.get_ligand(f"./NG1.mol2")
MG  = eh.get_metal('MG', charge=2)

stru.add(MG, chain_name='B', net_charge= 2, multiplicity=1) 
stru.add(ADP, chain_name='Y', net_charge= -3, multiplicity=1) 
stru.add(NG1, chain_name='Z', net_charge= -2, multiplicity=1) 

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")


eh.preparation.place_ligand(stru, 'Y', 1)
eh.preparation.place_ligand(stru, 'B', 1)
eh.preparation.place_ligand(stru, 'Z', 1, method='mole2', constraints=constraints)
sp.save_structure('assembled.pdb', stru)
exit( 0 ) 
#eh.preparation.place_ligand(stru, 'B', 1)
#
#sp.save_structure('temp.pdb', stru)
#stru = sp.get_structure('temp.pdb')

eh.dock_reactants( stru,
                    [stru.get('Z.1')],
                    constraints=constraints,
                    n_struct=20,
                    cluster_distance=2.5,
                    clash_cutoff=3,
                    cst_energy=1500,
                    )
