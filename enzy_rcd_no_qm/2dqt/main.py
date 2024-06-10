import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")
import enzy_htp as eh


eh.config['mole2.MONO'] = "/accre/arch/easybuild/software/Compiler/GCCcore/6.4.0/Mono/5.4.1.6/bin/mono"
eh.config['mole2.PROBE'] = 5.0
eh.config['mole2.INNER'] = 1.1
eh.config['mole2.IGNORE_HETATM'] = True


prefix='../../starting_points/2dqt/'



sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/filled.pdb')
CPD = lp.get_ligand(f"./CPD.mol2")

stru.add(CPD, chain_name='Z', net_charge=1, multiplicity=1) 

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")


eh.preparation.seed_ligand(stru, stru.get('Z.1'), method='mole2', constraints=constraints, minimize=True)
sp.save_structure('assembled.pdb', stru)

##
##sp.save_structure('temp.pdb', stru)
##stru = sp.get_structure('temp.pdb')
#
#eh.dock_reactants( stru,
#                    [stru.get('Z.1')],
#                    constraints=constraints,
#                    n_struct=20,
#                    cluster_distance=2.5,
#                    cst_energy=5000,
#                    box_size=100,
#                    max_sasa_ratio=1.0,
#                    clash_cutoff=50
#                    )
#
#sp.save_structure('final.pdb', stru)
