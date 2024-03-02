from openbabel import pybel
mol = next(pybel.readfile("mol2", "FAH.mol2"))


mol.calccharges()


mol.write("pdbqt", "FAH.pdbqt", overwrite=True)
