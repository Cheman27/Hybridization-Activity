import pandas as pd
from jupyter_jsmol import JsmolView
from ipywidgets import Layout, widgets, interact
from IPython.display import display, SVG
from IPython.display import Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from bs4 import BeautifulSoup
import math

# s orbital example code
def example_code_s_orbital():
    # set variable to add hydroges to our model. will call vairable later in another function. 
    eth_mol = Chem.AddHs(Chem.MolFromSmiles("C=C"))
    # takes the molecules connectivity to generate coordinates that can be used for 3D representations
    AllChem.EmbedMolecule(eth_mol)
    # converts smiles data to .mol file for better plotting later. 
    Chem.MolToMolFile(eth_mol, 'ethene.mol')
    # set variable to plot in Jsmol widget
    ec3Dso = JsmolView.from_file('ethene.mol', inline=True)
    # generate jsmol widget of our variable
    display(ec3Dso)
    # set background to white
    ec3Dso.script("background white")
    # refresh widget to prevent buggs
    ec3Dso.script("refresh")
    # set variable to make hybrid elements transparent to visually rich representations
    ec3Dso_transparent_line = ("select(atomno=1); lcaocartoon create 's'; color lcaocartoon transparent")
    # for loop to run transparent vairable twice. needed due to bug issues in python 
    for _ in range (2):
        ec3Dso.script(ec3Dso_transparent_line)


# p orbital
def example_code_p_orbital():
    eth_mol = Chem.AddHs(Chem.MolFromSmiles("C=C")) # set variable to add hydroges to our model. will call vairable later in another function.
    AllChem.EmbedMolecule(eth_mol) # takes the molecules connectivity to generate coordinates that can be used for 3D representations
    Chem.MolToMolFile(eth_mol, 'ethene.mol') # converts smiles data to .mol file for better plotting later. 
    
    ec3Dpo = JsmolView.from_file('ethene.mol', inline=True) # set variable to plot in Jsmol widget
    display(ec3Dpo) # generate jsmol widget of our variable
    ec3Dpo.script("background white") # set background to white
    ec3Dpo.script("refresh") # refresh widget to prevent buggs
    ec3Dpo_transparent_line = ("select(atomno=1); lcaocartoon create 'pz'; color lcaocartoon transparent") # set variable to make hybrid elements transparent to visually rich representations
    for _ in range(2): # for loop to run transparent vairable twice. needed due to bug issues in python 
        ec3Dpo.script(ec3Dpo_transparent_line)


# s and p orbital overlayed
def example_code_overlayed_orbital():
    eth_mol = Chem.AddHs(Chem.MolFromSmiles("C=C"))
    AllChem.EmbedMolecule(eth_mol)
    Chem.MolToMolFile(eth_mol, 'ethene.mol')
    ec3Doo = JsmolView.from_file('ethene.mol', inline=True)
    display(ec3Doo)
    ec3Doo.script("background white")
    ec3Doo.script("refresh")
    ec3Doo_transparent_line = ("select(atomno=1); lcaocartoon create 's'; lcaocartoon create 'pz'; color lcaocartoon transparent")
    for _ in range(2):
        ec3Doo.script(ec3Doo_transparent_line)


# an sp2 hybrid orbital
def example_code_hybridized_orbital():
    #add hydrogens from smiles information
    eth_mol = Chem.AddHs(Chem.MolFromSmiles("C=C"))
    #constrains distance within molecule
    AllChem.EmbedMolecule(eth_mol)
    #convert to .mol file
    Chem.MolToMolFile(eth_mol, 'ethene.mol')
    #set variable to view .mol file
    ec3Dho = JsmolView.from_file('ethene.mol', inline=True)
    #force display widget
    display(ec3Dho)
    #set background white for easier viewing
    ec3Dho.script("background white")
    #refresh the widget to avoid minor viewability issues
    ec3Dho.script("refresh")
    #set variable to generate one hybrid sp2 orbital
    ec3Dho_transparent_line = ("select(atomno=1); lcaocartoon create 'sp2c'; color lcaocartoon transparent")
    #make previous line run twice to workabout bug in function
    for _ in range(2):
        ec3Dho.script(ec3Dho_transparent_line)


#molecule 2D flat
def get_smiles(molecule_name):
    s_df = pd.read_csv("SMILES.txt", sep=",")
    s_df.columns = s_df.columns.str.strip()
    #return s_df
    molecule_smiles = s_df.loc[s_df['molecule'].str.lower() == molecule_name.lower(), 'SMILES'].values
    if molecule_smiles.size > 0:
        return molecule_smiles[0]
    else:
        print("Retry entry")
        return None

def view2D_flat(molecule_name):
    smiles = get_smiles(molecule_name)
    if smiles is not None:
        ma = Chem.MolFromSmiles(smiles)
        ma_with_h = Chem.AddHs(ma)
        return Draw.MolToImage(ma_with_h)


#formaldehyde 2D depth
def assign_visual_depth_multiple(mol, atom_idx, wedge_bonds=[], hash_bonds=[]):
    for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
        b_idx = bond.GetIdx()
        if b_idx in wedge_bonds:
            bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
        elif b_idx in hash_bonds:
            bond.SetBondDir(Chem.BondDir.BEGINDASH)
        else:
            bond.SetBondDir(Chem.BondDir.NONE)

def add_lone_pairs_oxygen_index(svg, oxygen_idx=1):
    soup = BeautifulSoup(svg, "xml")

    # Find the group for the oxygen atom (usually atom-1)
    oxygen_group = soup.find("g", {"class": f"atom-{oxygen_idx}"})
    if not oxygen_group:
        print("Could not find oxygen group.")
        return svg

    # Try to find the 'circle' or 'text' to extract coordinates
    circle = oxygen_group.find("circle")
    if circle:
        x = float(circle["cx"])
        y = float(circle["cy"])
    else:
        text = oxygen_group.find("text")
        if text:
            x = float(text["x"])
            y = float(text["y"])
        else:
            print("Could not find position for oxygen.")
            return svg

    # Add lone pair dots above
    dxs = [-3, 3]
    dy = 6
    dot_radius = 1.2

    for dx in dxs:
        dot = soup.new_tag("circle", cx=str(x + dx), cy=str(y - dy),
                           r=str(dot_radius), fill="red")
        oxygen_group.append(dot)

    return str(soup)

# Final combined function
def view2D_depth_formaldehyde():
    mol = Chem.AddHs(Chem.MolFromSmiles("C=O"))
    AllChem.EmbedMolecule(mol)
    AllChem.Compute2DCoords(mol)

    idx_C = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C'][0]
    idx_O = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'O'][0]
    C_bonds = [b.GetIdx() for b in mol.GetAtomWithIdx(idx_C).GetBonds()]
    assign_visual_depth_multiple(mol, idx_C, wedge_bonds=[C_bonds[0]], hash_bonds=[C_bonds[0]])

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    
    # Get drawing coordinates of oxygen (in SVG space)
    draw_coords = drawer.GetDrawCoords(idx_O)
    svg_x = draw_coords.x
    svg_y = draw_coords.y

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    soup = BeautifulSoup(svg, "xml")

    # Add lone pair dots
    dxs = [-6, 6]
    dy = 25
    dot_radius = 3

    # Lone pair 1, above oxygen
    for dx in [-6, 6]:
        dot = soup.new_tag("circle",
                           cx=str(svg_x + dx + 3),
                           cy=str(svg_y - 24),
                           r="2", fill="red")
        soup.svg.append(dot)
    
    # Lone pair 2, below oxygen
    for dx, dy in [(9, -2), (12, -13)]:
        dot = soup.new_tag("circle",
                           cx=str(svg_x + dx + 10),
                           cy=str(svg_y + dy + 17),
                           r="2", fill="red")
        soup.svg.append(dot)

    display(SVG(str(soup)))


#acetonitrile 2D depth
def add_lone_pairs_nitrogen_index(svg, oxygen_idx=1):
    soup = BeautifulSoup(svg, "xml")

    # Find the group for the oxygen atom (usually atom-1)
    nitrogen_group = soup.find("g", {"class": f"atom-{nitrogen_idx}"})
    if not nitrogen_group:
        print("Could not find nitrogen group.")
        return svg

    # Try to find the 'circle' or 'text' to extract coordinates
    circle = nitrogen_group.find("circle")
    if circle:
        x = float(circle["cx"])
        y = float(circle["cy"])
    else:
        text = nitrogen_group.find("text")
        if text:
            x = float(text["x"])
            y = float(text["y"])
        else:
            print("Could not find position for nitrogen.")
            return svg

    # Add lone pair dots above
    dxs = [-3, 3]
    dy = 6
    dot_radius = 1.2

    for dx in dxs:
        dot = soup.new_tag("circle", cx=str(x + dx), cy=str(y - dy),
                           r=str(dot_radius), fill="red")
        nitrogen_group.append(dot)

    return str(soup)

# Final combined function
def view2D_depth_acetonitrile():
    mol = Chem.AddHs(Chem.MolFromSmiles("CC#N"))
    AllChem.EmbedMolecule(mol)
    AllChem.Compute2DCoords(mol)
    
    # Get atom indices
    idx_C = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C'][0]
    idx_C2 = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C'][0]
    idx_N = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'N'][0]

    # Bond indices around carbon1
    C_bonds = [b.GetIdx() for b in mol.GetAtomWithIdx(idx_C).GetBonds()]
    assign_visual_depth_multiple(mol, idx_C, wedge_bonds=[C_bonds[0]], hash_bonds=[C_bonds[0]])

    #Bond indicies around carbon2
    C2_bonds = [b.GetIdx() for b in mol.GetAtomWithIdx(idx_C2).GetBonds()]
    assign_visual_depth_multiple(mol, idx_C, wedge_bonds=[C_bonds[1]], hash_bonds=[C_bonds[3]])

    # Bond indices around nitrogen
    N_bonds = [b.GetIdx() for b in mol.GetAtomWithIdx(idx_N).GetBonds()]
    
    # Keep N–C bond plain. Assume N_bonds[0] is N–C, others are N–H
    assign_visual_depth_multiple(mol, idx_N, wedge_bonds=[], hash_bonds=[N_bonds[0], N_bonds[0]])

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    
    # Get drawing coordinates of nitrogen (in SVG space)
    draw_coords = drawer.GetDrawCoords(idx_N)
    svg_x = draw_coords.x
    svg_y = draw_coords.y

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    soup = BeautifulSoup(svg, "xml")

    # Add lone pair dots
    dxs = 0
    dy = [-6, 6]
    dot_radius = 3

    # Lone pair 1, adjacent nitrogen
    for dy in [-6, 6]:
        dot = soup.new_tag("circle",
                           cx=str(svg_x - 15),
                           cy=str(svg_y + dy),
                           r="2", fill="blue")
        soup.svg.append(dot)
    display(SVG(str(soup)))


#pyrrole 2D depth
def add_lone_pairs_nitrogen_index(svg, oxygen_idx=1):
    soup = BeautifulSoup(svg, "xml")

    # Find the group for the oxygen atom (usually atom-1)
    nitrogen_group = soup.find("g", {"class": f"atom-{nitrogen_idx}"})
    if not nitrogen_group:
        print("Could not find nitrogen group.")
        return svg

    # Try to find the 'circle' or 'text' to extract coordinates
    circle = nitrogen_group.find("circle")
    if circle:
        x = float(circle["cx"])
        y = float(circle["cy"])
    else:
        text = nitrogen_group.find("text")
        if text:
            x = float(text["x"])
            y = float(text["y"])
        else:
            print("Could not find position for nitrogen.")
            return svg


# Final combined function
def view2D_depth_pyrrole():
    mol = Chem.AddHs(Chem.MolFromSmiles("C1=CNC=C1"))
    AllChem.EmbedMolecule(mol)
    AllChem.Compute2DCoords(mol)

    idx_N = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'N'][0]
    N_bonds = [b.GetIdx() for b in mol.GetAtomWithIdx(idx_N).GetBonds()]
    
    # Keep N–C bond plain. Assume N_bonds[0] is N–C, others are N–H
    assign_visual_depth_multiple(mol, idx_N, wedge_bonds=[], hash_bonds=[])

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    
    # Get drawing coordinates of nitrogen (in SVG space)
    draw_coords = drawer.GetDrawCoords(idx_N)
    svg_x = draw_coords.x
    svg_y = draw_coords.y

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    soup = BeautifulSoup(svg, "xml")

    # Add lone pair dots
    dxs = 0
    dy = [-6, 6]
    dot_radius = 3

    # Lone pair 1, adjacent nitrogen
    for dx, dy in [(-6, 6), (1, -6)]:
        dot = soup.new_tag("circle",
                           cx=str(svg_x + dx + 20),
                           cy=str(svg_y + dy + 10),
                           r="2", fill="blue")
        soup.svg.append(dot)
    display(SVG(str(soup)))


#methylamine 3D static
def view3D_static_methylamine():
    display(Image(filename='3D_static_methylamine.png'))


#formaldehyde 3D static
def view3D_static_formaldehyde():
    display(Image(filename='3D_static_formaldehyde.png'))


#acetonitrile 3D static
def view3D_static_acetonitrile():
    display(Image(filename='3D_static_acetonitrile.png'))


#pyrrole 3D static
def view3D_static_pyrrole():
    display(Image(filename='3D_static_pyrrole.png'))


#methylamine 3D interactive 
def view3D_int_methylamine(): 
    macords3D = JsmolView.from_file('methylamine.xyz', inline=True)
    #need to find command to force view 
    display(macords3D)
    macords3D.script("background white")
    macords3D.script("refresh")
    macords3D.script("select(atomno=1); lcaocartoon lonepair 'lpa'")
    transparent_line = "select(atomno=1); lcaocartoon create 'sp3a'; lcaocartoon create 'sp3b'; lcaocartoon create 'sp3c'; lcaocartoon create 'lp'; color lcaocartoon translucent"
    for _ in range (2):
        macords3D.script(transparent_line)


#formaldehyde 3D interactive
def view3D_int_formaldehyde():
    fmol = Chem.AddHs(Chem.MolFromSmiles("C=O"))
    AllChem.EmbedMolecule(fmol)
    Chem.MolToMolFile(fmol, "formaldehyde.mol")

    fcords3D = JsmolView.from_file('formaldehyde.mol', inline=True)
    #force view 
    display(fcords3D)
    fcords3D.script("background white")
    fcords3D.script("refresh")
    #need to add lonepairs at sp2 locations
    fcords3D.script("select(atomno=2); lcaocartoon lonepair 'sp2b'; lcaocartoon lonepair 'sp2c'")
    fcords3D.script("select(atomno=2); lcaocartoon create 'pz'")
    fcords3D.script("select(atomno=1); lcaocartoon create 'pz'")
    #need to make next line run twice, to make lobes transparent
    transparent_f = "select(atomno=2); lcaocartoon create 'sp2a'; lcaocartoon create 'sp2b'; lcaocartoon create 'sp2c'; color lcaocartoon translucent"
    for _ in range (2):
        fcords3D.script(transparent_f)


#acetonitrile 3D interactive
def view3D_int_acetonitrile():
    anmol = Chem.AddHs(Chem.MolFromSmiles("CC#N"))
    AllChem.EmbedMolecule(anmol)
    Chem.MolToMolFile(anmol, "acetonitrile.mol")

    ancords3D = JsmolView.from_file('acetonitrile.mol', inline=True)
    #force view 
    display(ancords3D)
    ancords3D.script("background white")
    ancords3D.script("refresh")
    #need to add lonepairs at sp location
    ancords3D.script("select(atomno=3); lcaocartoon lonepair 'spb'")
    #limitation: dificult to show the sigma bond hybridized, and the lone pair hybridized 
    transparent_n = "select(atomno=3); lcaocartoon create 'px'; lcaocartoon create 'py'; lcaocartoon create 'spa'; lcaocartoon create 'spb'; color lcaocartoon translucent"
    transparent_c = "select(atomno=2); lcaocartoon create 'px'; lcaocartoon create 'py'; color lcaocartoon translucent"
    for _ in range (2):
         ancords3D.script(transparent_n)
    for _ in range (2):
        ancords3D.script(transparent_c)


#acetonitrile 3D interactive
def view3D_int_acetonitrile():
    anmol = Chem.AddHs(Chem.MolFromSmiles("CC#N"))
    AllChem.EmbedMolecule(anmol)
    Chem.MolToMolFile(anmol, "acetonitrile.mol")

    ancords3D = JsmolView.from_file('acetonitrile.mol', inline=True)
    #force view 
    display(ancords3D)
    ancords3D.script("background white")
    ancords3D.script("refresh")
    #need to add lonepairs at sp location
    ancords3D.script("select(atomno=3); lcaocartoon lonepair 'spb'")
    #limitation: dificult to show the sigma bond hybridized, and the lone pair hybridized 
    transparent_n = "select(atomno=3); lcaocartoon create 'px'; lcaocartoon create 'py'; lcaocartoon create 'spa'; lcaocartoon create 'spb'; color lcaocartoon translucent"
    transparent_c = "select(atomno=2); lcaocartoon create 'px'; lcaocartoon create 'py'; color lcaocartoon translucent"
    for _ in range (2):
         ancords3D.script(transparent_n)
    for _ in range (2):
        ancords3D.script(transparent_c)


#pyrrole 3D interactive
def view3D_int_pyrrole():
    pmol = Chem.AddHs(Chem.MolFromSmiles("C1=CNC=C1"))
    AllChem.EmbedMolecule(pmol)
    Chem.MolToMolFile(pmol, "pyrrole.mol")

    pcords3D = JsmolView.from_file('pyrrole.mol', inline=True)
    #force view 
    display(pcords3D)
    pcords3D.script("background white")
    pcords3D.script("refresh")
    p_translucent = ("select(atomno=3); lcaocartoon create 'pz'; color lcaocartoon translucent")
    for _ in range (2):
        pcords3D.script(p_translucent)