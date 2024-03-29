import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh

eh.config['mole2.MONO'] = "/accre/arch/easybuild/software/Compiler/GCCcore/6.4.0/Mono/5.4.1.6/bin/mono"
eh.config['mole2.PROBE'] = 5.0
eh.config['mole2.INNER'] = 1.05
eh.config['mole2.IGNORE_HETATM'] = True


prefix='../../starting_points/1nbh/'

sp = eh.PDBParser()
lp = eh.Mol2Parser()

stru = sp.get_structure(f'{prefix}/full_chain_renumbered.pdb')

SAM = lp.get_ligand(f"{prefix}/SAM.mol2")
ACT = lp.get_ligand(f"{prefix}/ACT.mol2")

stru.add(SAM, chain_name="Y", net_charge= 1, multiplicity=1)
stru.add(ACT, chain_name='Z', net_charge=-1, multiplicity=1)

constraints = eh.structure.structure_constraints_from_xml(stru, "constraints.xml")

eh.preparation.place_ligand(stru, 'Y', 1)
eh.preparation.place_ligand(stru, 'Z', 1, method='mole2', constraints=constraints)

sp.save_structure('assembled.pdb', stru)


eh.dock_reactants( stru,
                    [stru.get('Z.1')],
                    constraints=constraints,
                    n_struct=100,
                    cluster_distance=2.5,
                    clash_cutoff=3,
                    cst_energy=2000,
                    )


sp.save_structure('result.pdb', stru)
