import os
import shutil
from pymol import cmd, stored
from pathlib import Path
from openbabel import pybel


ligands = ["ADP.mol2",  "Mg.mol2", "NG1.mol2"]


code = Path('.').parent.absolute().stem


vina="/Users/chrisjurich/projects/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina"

cmd.load(f'../../starting_points/{code}/full_chain_renumbered.pdb')

stored.holder = list()

cmd.iterate('all', 'stored.holder.append(resi)')
curr = max(map(int, stored.holder)) + 1


cmd.save('start.pdb')
os.system('pdb2pqr30 --keep-chain start.pdb start.pdb')

shutil.copy('start.pdb', 'intermediate.pdb')

docked = list()

for lidx, ll in enumerate(ligands):

    os.system('obabel -i pdb intermediate.pdb -o pdbqt -O intermediate.pdbqt -p 7 -xr 1> /dev/null ')
    
    mol = next(pybel.readfile("mol2", ll))
    mol.calccharges()
    print(ll)
    lpath = Path(ll)
    print(lpath.with_suffix('.pdbqt'))
    ll = str(lpath.with_suffix('.pdbqt'))
    mol.write("pdbqt", ll, overwrite=True)

    ligands[lidx] = ll

    config_txt = f"config{lidx}.txt"
    resn = Path(ll).stem

    cmd.delete('all')
    cmd.fetch(code)
    cmd.remove('solvent')
    [com_x, com_y, com_z] = cmd.centerofmass(f"resn {resn}")

    fh = open(config_txt, 'w')

    fh.write(f"""center_x  = {com_x} 
center_y  = {com_y} 
center_z  = {com_z} 
size_x    = 10.0
size_y    = 10.0
size_z    = 10.0""")

    fh.close()

    cmd_str:str=f"{vina} --exhaustiveness 40 --receptor intermediate.pdbqt --config {config_txt} --ligand {ll}"
    print(cmd_str)
    os.system(cmd_str)
    tpath = Path(ll)
    outfile = f"{resn}_out.pdbqt"

    cmd.delete('all')
    cmd.load(outfile)
    cmd.split_states(f"{resn}_out")
    
    for on in cmd.get_object_list()[1:]:
        cmd.delete(on)

    cmd.alter('all', 'chain=""')
    cmd.alter('all', f'resi="{curr}"')
    curr += 1
    lig_stru = f"{resn}_out.pdb"
    
    cmd.save(lig_stru)
    
    docked.append(lig_stru)
    cmd.delete('all')
    
   
    cmd.load('start.pdb')
    for ls in docked:
        cmd.load(ls)
    
    cmd.save('intermediate.pdb')
    
