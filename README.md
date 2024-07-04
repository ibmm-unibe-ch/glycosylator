
![logo](docs/source/_static/_resources/logo_large.png)

# A Python framework for the rapid modeling of glycans

Built using [BuildAMol](https://github.com/NoahHenrikKleinschmidt/buildamol) Glycosylator is a Python framework for the identification, modeling and
modification of glycans. Glycosylator can build atomic models of glycans and can glycosylate proteins and membranes. Glycosylator can perform conformational optimization to minimize clashes in glycosylated structures or to sample alternative conformations for individual glycans. Glycosylator supports a variety of file-formats and can work hand in hand with other libraries such as RDKit to faciliate research workflows. 


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

<!-- Please cite:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3097-6 -->


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
