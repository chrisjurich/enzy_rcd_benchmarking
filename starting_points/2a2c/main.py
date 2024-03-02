import sys
sys.path.insert(0, "/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/")

import enzy_htp as eh



session = eh.interface.pymol.new_session()

eh.interface.pymol.general_cmd(session, [

('fetch', 'NG1'),
('save', 'NG1.mol2'),
('delete', 'all'),
('fetch', 'ADP'),
('save', 'ADP.mol2'),
])

lp = eh.Mol2Parser()



eh.interface.moe.protonate('NG1.mol2')
eh.interface.moe.protonate('ADP.mol2')

