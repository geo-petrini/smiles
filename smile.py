from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor


mols = []
mols.append( Chem.MolFromSmiles('Cc1ccccc1')  )
mols.append( Chem.MolFromSmiles('CC(=O)OCC[N+](C)(C)C') )
mols.append( Chem.MolFromSmiles('C1=CC=C(C=C1)C=O') )
mols.append( Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O') )
mols.append( Chem.AddHs(Chem.MolFromSmiles('O') )) #acqua
mols.append( Chem.AddHs(Chem.MolFromSmiles('C') ) )#metano
mols.append( Chem.AddHs(Chem.MolFromSmiles('C=C') )) #ethano
mols.append( Chem.AddHs( Chem.MolFromSmiles('O=CCCO')) ) #ethano
mols.append( Chem.MolFromSmiles('[H]O[H]') ) #acqua

# m = Chem.MolFromSmiles('C(=C(Cl)Cl)(Cl)Cl ') #terracloro


# m = Chem.MolFromSmiles('C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O', params)

for i, m in enumerate(mols):
    m.Debug()

    rdDepictor.Compute2DCoords(m)
    rdDepictor.StraightenDepiction(m)

    # img = Draw.MolToImage(m)
    # im1 = img.save("mol.png")
    # Draw.MolToFile(m, "mol.png")
    # for atom in m.GetAtoms():
    #     atom.SetAtomMapNum(atom.GetIdx())

    d = Draw.MolDraw2DCairo(300, 300)
    # d = Draw.MolDraw2DSVG(300, 300)
    dopts = d.drawOptions()
    dopts.addAtomIndices  = True
    dopts.explicitMethyl = True
    dopts.includeRadicals = False
    dopts.useComplexQueryAtomSymbols = False
    legend = f'{Chem.MolToSmiles(m)}'
    # legend += f'\naddAtomIndices: {dopts.addAtomIndices}\nexplicitMethyl: {dopts.explicitMethyl}'
    d.DrawMolecule(m, legend=legend)
    # d.FinishDrawing()
    d.WriteDrawingText(f"mol_{i}.png")
    # with open(f"mol_{i}.svg", 'w') as f:
    #     f.write( d.GetDrawingText() )




