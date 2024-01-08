import os
import sys
import shutil
import enzy_htp as eh
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from collections import defaultdict
sys.path.insert(0, '/Library/modeller-10.4//modlib/')
from modeller import *
from modeller.automodel import *    # Load the AutoModel class
from pymol import cmd
from string import ascii_uppercase 

sys.setrecursionlimit(10000)
def renumber(session, selection='all', start=1, startsele=None, quiet=1):
    '''
DESCRIPTION

    Set residue numbering (resi) based on connectivity.

ARGUMENTS

    selection = string: atom selection to renumber {default: all}

    start = integer: counting start {default: 1}

    startsele = string: residue to start counting from {default: first in
    selection}
    '''
    start, quiet = int(start), int(quiet)
    model = session.cmd.get_model(selection)
    session.cmd.iterate(selection, 'next(atom_it).model = model',
            space={'atom_it': iter(model.atom), 'next': next})
    if startsele is not None:
        startidx = session.cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    session.cmd.alter(selection, 'resi = next(atom_it).resi',
            space={'atom_it': iter(model.atom), 'next': next})
    if not quiet:
        print(' Renumber: range (%d to %d)' % tuple(minmax))


parser = eh.PDBParser()

structure = parser.get_structure("./raw_pdb.pdb")

eh.preparation.remove_solvent(structure)

pp_chains = list()
for cc in structure.chains:
    if cc.is_polypeptide():
        pp_chains.append( cc.name )

parser.save_structure("./raw_pdb_no_wat.pdb", structure)

pp_sele = " or ".join(map(lambda cc: f"( chain {cc})", pp_chains))
session = eh.interface.pymol.new_session()
eh.interface.pymol.general_cmd(session, [
    ('load', "./raw_pdb_no_wat.pdb"),
    ("remove", f"not ({pp_sele})"),
    ("save", "./raw_pdb_no_wat.pdb")
])
df = eh.interface.pymol.collect(session, "memory", "chain resi resn".split())

raw_residues = set()

for i, row in df.iterrows():
    raw_residues.add(
        (row['chain'],int(row['resi']), row['resn'], False)
    )

parser = PDBParser()

structure = parser.get_structure("stru", "./raw_pdb.pdb")

for dd in structure.header['missing_residues']:
    raw_residues.add(
        (dd['chain'], dd['ssseq'], dd['res_name'], True )
    )

grouped_residues=defaultdict(list)    

for rr in raw_residues:
    grouped_residues[rr[0]].append(
        {'resi': rr[1], 'resn':rr[2], 'missing':rr[3]}
    )


curr_seq = ''
target_seq = ''

res_keys = list()
for k in sorted(list(grouped_residues.keys())):
    v = grouped_residues[k]
    v.sort(key=lambda row: row['resi'])
    grouped_residues[k] = v
    temp = ''
    for vv in v:
        if vv['resn'] in 'ADP NG1 MG NA'.split():
            curr_seq += '?'
            target_seq += '?'
        else:
            one_letter = eh.chemical.residue.convert_to_one_letter(vv['resn'])
            target_seq += one_letter
            if vv['missing']:
                curr_seq += '-'
            else:
                
                res_keys.append((k, vv['resi']))
                curr_seq += one_letter
        
    curr_seq += "/"
    target_seq += "/"


curr_seq = curr_seq[:-1]
target_seq = target_seq[:-1]

curr_seq += '*'
target_seq += '*'

offset_mapper = dict()

for k, v in grouped_residues.items():
    
    for ridx, r_dict in enumerate(v):
        offset_mapper[k] = r_dict['resi']
        break

if curr_seq != target_seq:
    fh = open('alignment.ali', 'w')
    
    fh.write(
    f""">P1;raw_pdb_no_wat
structureX:raw_pdb_no_wat:{res_keys[0][1]}:{res_keys[0][0]}:{res_keys[-1][1]}:{res_keys[-1][0]}:undefined:undefined:-1.00:-1.00
{curr_seq}
>P1;raw_pdb_no_wat_fill
sequence:::::::::
{target_seq}"""
    )
    
    fh.close()
    
    log.none()
    env = Environ()
    
    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']
    
    a = LoopModel(env, alnfile = 'alignment.ali',
                  knowns = 'raw_pdb_no_wat', sequence = 'raw_pdb_no_wat_fill')
    a.starting_model= 1
    a.ending_model  = 1
    
    a.loop.starting_model = 1
    a.loop.ending_model   = 2
    a.loop.md_level       = refine.fast
    
    a.make()
    
    
    shutil.move("raw_pdb_no_wat_fill.B99990001.pdb", "full_chain.pdb")

    for tk in """raw_pdb_no_wat_fill.BL00010001.pdb
    2a2c.cif
    raw_pdb_no_wat_fill.B99990001.pdb
    raw_pdb_no_wat_fill.BL00020001.pdb
    raw_pdb_no_wat_fill.D00000001
    raw_pdb_no_wat_fill.DL00010001
    raw_pdb_no_wat_fill.DL00020001
    raw_pdb_no_wat_fill.IL00000001.pdb
    raw_pdb_no_wat_fill.ini
    raw_pdb_no_wat_fill.lrsr
    raw_pdb_no_wat_fill.rsr
    raw_pdb_no_wat_fill.sch
    raw_pdb_no_wat_fill.V99990001""".split():
        if os.path.exists(tk):
            os.remove(tk)
    
    session = eh.interface.pymol.new_session()
    
    eh.interface.pymol.general_cmd(session, [
        ('load', 'full_chain.pdb'),
        ('load', 'raw_pdb_no_wat.pdb'), 
        ('align', 'full_chain', 'raw_pdb_no_wat'),
        ('delete', 'raw_pdb_no_wat'),
    ])

    chain_names = sorted(list(offset_mapper.keys()))
    
    chain_remapper = dict()

    for cnidx, cn in enumerate(chain_names):
        renumber(session, selection=f'chain {ascii_uppercase[cnidx]}', start=offset_mapper[cn])
           
        if ascii_uppercase[cnidx] != cn:
            chain_remapper[ascii_uppercase[cnidx]] = cn            

    if chain_remapper:
        session.cmd.alter('all', 'chain=c_mapper[chain]', space={'c_mapper':chain_remapper})            

    
    session.cmd.save('full_chain_renumbered.pdb')

else:
    shutil.copy("raw_pdb_no_wat.pdb", "full_chain_renumbered.pdb")

