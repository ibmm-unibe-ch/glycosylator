import base as base
import glycosylator as gl


def test_shield1():
    protein = gl.Protein.load(
        base.HOME + "/docs/source/examples/files/protein_optimized.pkl"
    )
    for i, glycan in enumerate(protein.glycans):
        glycan.id = f"glycan-{i}"

    # setup the shield simulator
    shield = gl.GlycoShield(protein)

    # now run the simulation
    # (and visualize the conformations)
    out = shield.simulate(
        angle_step=50, save_conformations_to=None, visualize_conformations=False
    )
    pass
