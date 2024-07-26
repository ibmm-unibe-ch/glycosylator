
![logo](docs/source/_static/_resources/logo_large.png)

# A Python framework for the rapid modeling of glycans
[![Made with Python](https://img.shields.io/badge/Python->=3.10-blue?logo=python&logoColor=white)](https://python.org "Go to Python homepage")
[![code style - black](https://img.shields.io/badge/code_style-black-black)](https://black.readthedocs.io/ "Go to Black homepage")
[![Documentation Status](https://readthedocs.org/projects/glycosylator/badge/?version=latest)](https://glycosylator.readthedocs.io/en/latest/?badge=latest)
[![Check out - Tutorials](https://img.shields.io/badge/check_out-Tutorials-e61882)](https://glycosylator.readthedocs.io/en/latest/tutorials.html)
[![PyPI version](https://badge.fury.io/py/glycosylator.svg)](https://badge.fury.io/py/glycosylator)

Built using [BuildAMol](https://github.com/NoahHenrikKleinschmidt/buildamol) Glycosylator is a Python framework for the identification, modeling and
modification of glycans. Glycosylator can build atomic models of glycans and can glycosylate proteins and membranes. Glycosylator can perform conformational optimization to minimize clashes in glycosylated structures or to sample alternative conformations for individual glycans. Glycosylator supports a variety of file-formats and can work hand in hand with other libraries such as RDKit to faciliate research workflows. 

Running the GUI
---------------

[![Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://glycosylator.streamlit.app)

Glycosylator provides a graphical user interface in the form of a Streamlit Web Application. This app supports modeling new glycans and glycosylating proteins and memrbanes as well as modeling the shielding effect. Of course, the library itself is more versatile and some things may not be doable via the GUI. The app can be reached via [glycosylator.streamlit.app](glycosylator.streamlit.app) or via the _Open in Streamlit_ badge at the top. If you would like to run the app locally you clone [the app's Github Repository](https://github.com/NoahHenrikKleinschmidt/Streamlit-Glycosylator). 

Installing Glycosylator
-----------------------
Glycosylator is distributed via the Python Package Index and can be installed via:

```bash
pip install glycosylator
```


What can Glycosylator do?
-------------------------
Here's a list of things that you could do with Glycosylator: 

- Model new glycans
    - From IUPAC strings
    - From GlyCosmos IDs
    - Manually through fragment-based assembly
    - Create non-standard sugars by adding functional groups or other small molecule fragments

- Modify existing glycans
  - Extend small glycans to larger ones
  - Trim large glycans to smaller ones
  - Modify individual sugar residues
  - Sample or optimize glycan conformations
  - Rename atoms/residues to conform to different forcefield conventions

- Obtain structural data on glycans
  - Get data for angles / dihedrals
  - Get compositional data for sugar residues

- Glycosylate Biomolecules
  - Glycosylate Proteins at specific residues or search by Sequons
  - Glycosylate Membranes at specific lipids such as Ceramide
  - Glycosylate arbitrary molecular scaffolds in pretty much any way

- Sample or optimize glycan conformations on the surface of scaffold molecules
- Simulate the shielding effect of glycans on scaffold molecules


Building a glycan from IUPAC
----------------------------

Generating an atomic model from a IUPAC glycan string can be as easy as:

```python
import glycosylator as gl

# this is GlyCosmos entry G09568JT
iupac = "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"

glycan = gl.glycan(iupac)
glycan.to_pdb("my_glycan.pdb")

glycan.show3d()
```

> ![](docs/source/_static/_resources/glycan_light.gif)
> Note: This visualization was hand-made in ChimeraX from the PDB file, so the actual output of show3d will not be quite this pretty, but it will be interactive... 

Glycosylating a protein
-----------------------
To glycosylate a protein we can use a simple snippet like the one below. Check out the tutorials in the documentation for more info on glycosylation and how to optimize a glycosylated protein afterward.

```python
import glycosylator as gl

# get some protein scaffold
protein = gl.Protein.from_pdb("my_protein.pdb")

# get some glycan (maybe we already have one as a PDB file)
glycan = gl.Glycan.from_pdb("my_glycan.pdb")

# glycosylate the protein at N-glycosylation sites
glycoprotein = gl.glycosylate(protein, glycan, sequon="N-linked")

glycoprotein.to_pdb("my_glycoprotein.pdb")
```

> ![](docs/source/_static/_resources/prot_glyco.gif)

Simulating Glycan Shielding
---------------------------
To simulate the shielding effect of glycans on the surface of proteins we can use the `quickshield` function or the `GlycoShield` class. These will perform torsional sampling around random subsets of edges in order to obtain a great number of conformations in a relatively short time. Here is how: 

```python
# using a previously glycosylated protein
protein = gl.Protein.load("my_glycoprotein.pkl")

# perform a quick simulation of glycan shielding
# (with very rough settings)
shield = gl.quickshield(protein, angle_step=100, repeats=1, save_conformations_to="./glycoshield_simulation")

```

> ![](docs/source/_static/_resources/glycoshield.gif)


<!-- Please cite:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3097-6 -->


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
