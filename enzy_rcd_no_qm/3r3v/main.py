import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh


eh.config['mole2.MONO'] = "/accre/arch/easybuild/software/Compiler/GCCcore/6.4.0/Mono/5.4.1.6/bin/mono"
eh.config['mole2.PROBE'] = 5.0
eh.config['mole2.INNER'] = 1.05
eh.config['mole2.IGNORE_HETATM'] = True


prefix='../../starting_points/3r3v/'



sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/filled.pdb')
FAH = lp.get_ligand(f"{prefix}/FAH.mol2")

stru.add(FAH, chain_name='Z', net_charge= -1, multiplicity=1) 

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")


eh.preparation.seed_with_constraints(stru.get('Z.1'), constraints)
sp.save_structure('assembled.pdb', stru)

eh.dock_reactants( stru,
                    [stru.get('Z.1')],
                    constraints=constraints,
                    grid_width=100,
                    box_size=100,
                    n_struct=20,
                    transform_temperature=30,
                    cluster_distance=2.5,
                    clash_cutoff=3,
                    cst_energy=1500,
                    )

sp.save_structure('3r3v_enzy_rcd_no_qm.pdb', stru)
