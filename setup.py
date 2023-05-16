import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="glycosylator",
    version="2.0.38-alpha",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A python package for structural simulation of protein glycosylation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={
        "glycosylator.resources": ["*.pkl"],
    },
    url="TBA",
    packages=[
        "glycosylator",
        "glycosylator.core",
        "glycosylator.resources",
        "glycosylator.utils",
        "glycosylator.structural",
        "glycosylator.graphs",
        "glycosylator.optimizers",
    ],
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "networkx",
        "biopython",
        "pdbecif",
        "periodictable",
        "plotly",
        "alive_progress",
        "pyswarms",
        # "openbabel",
        # "mdtraj",
        # "prody",
        "scikit-learn",
        "nglview",
        "ipywidgets<8.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
)
