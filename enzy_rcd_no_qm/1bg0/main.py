import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh

prefix='../../starting_points/1bg0/'

eh.config['mole2.MONO'] = "/accre/arch/easybuild/software/Compiler/GCCcore/6.4.0/Mono/5.4.1.6/bin/mono"
eh.config['mole2.PROBE'] = 5.0
eh.config['mole2.INNER'] = 1.05
eh.config['mole2.IGNORE_HETATM'] = False

#TODO(CJ): need to put the other place_ligand() method

sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/full_chain_renumbered.pdb')

ADP = lp.get_ligand(f"./ADP.mol2")
DAR = lp.get_ligand(f"./DAR.mol2")
NO3 = lp.get_ligand(f"./NO3.mol2")
MG  = eh.get_metal('MG', charge=2)

stru.add(MG , chain_name="B", net_charge= 2, multiplicity=1)
stru.add(ADP, chain_name='X', net_charge=-3, multiplicity=1)
stru.add(NO3, chain_name='Y', net_charge=-1, multiplicity=1) 
stru.add(DAR, chain_name='Z', net_charge= 1, multiplicity=1) 

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")

eh.preparation.seed_with_transplants(stru.get('X.1'), 'atom_names')
eh.preparation.seed_with_transplants(stru.get('Z.1'), 'mcs', similarity_cutoff=0.20)
eh.preparation.seed_with_constraints(stru.get('Y.1'), [constraints[-1]])
eh.preparation.seed_with_transplants(stru.get('B.1'), 'atom_names')

sp.save_structure('temp.pdb', stru)


eh.dock_reactants( stru,
                    [stru.get('Y.1'), stru.get('Z.1')],
                    constraints=constraints,
                    n_struct=10,
                    cluster_distance=2.5,
                    clash_cutoff=3,
                    cst_energy=1500,
                    )

sp.save_structure('1bg0_enzy_rcd_no_qm.pdb', stru)
