import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="glycosylator",
    version="4.4.0",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A python package for structural simulation of protein glycosylation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={
        "glycosylator.resources": ["*.pkl", "*.json"],
        "glycosylator.resources.icons": ["*.png"],
    },
    url="https://github.com/ibmm-unibe-ch/glycosylator/",
    packages=[
        "glycosylator",
        "glycosylator.core",
        "glycosylator.utils",
        "glycosylator.resources",
        "glycosylator.optimizers",
        "glycosylator.resources.icons",
        "glycosylator.resources.names",
    ],
    install_requires=[
        "biobuild",
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
        "gym",
        "pubchempy",
        "tabulate",
        "scikit-learn",
        "ipywidgets<8.0.1",
    ],
    optional_requires={
        "openbabel": ["openbabel"],
        "rdkit": ["rdkit"],
        "glycowork": ["glycowork"],
        "full": ["openbabel", "rdkit", "glycowork"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
)
