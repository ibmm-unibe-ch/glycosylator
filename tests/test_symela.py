import glycosylator as gl


def test_symelas_protein():
    DEFAULT_IUPAC = (
        "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)"
        "[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]"
        "Man(b1-4)GlcNAc(b1-4)GlcNAc"
    )
    glycan = gl.glycan(DEFAULT_IUPAC)

    protein = gl.Protein.from_pdb(
        "/Users/noahhk/Downloads/fold_mev_nag_model_0_protein.pdb"
    )
    protein.reindex()

    n_sites = protein.find_n_linked_sites()
    sites = []
    for chain, residues in n_sites.items():
        protein.apply_standard_bonds_for(*residues)
        for residue in residues:
            nd2 = residue.get_atom("ND2")
            protein.add_hydrogens(nd2)
            for h in nd2.get_hydrogens():
                protein.rename_atom(h, f"HD{h.name[-2:]}")
            sites.append(residue)

    print("Starting glycosylation...")
    protein.attach(
        glycan,
        link=gl.get_linkage("ASN-glyco"),
        residues=sites,
        chain="new",
    )

    out = gl.optimize_glycans(protein, seed=0, maxiter=10)

    # print("Glycan attached, starting scaffold graph generation...")
    # graph, edges = gl.make_scaffold_graph(
    #     protein,
    #     only_clashing_glycans=False,
    #     include_root=True,
    #     include_n_ancestors=2,
    # )

    # print("Starting optimization environment setup...")
    # env = gl.optimizers.DistanceRotatron(graph, edges)
    # env = gl.optimizers.ScaffoldRotatron(env)

    # sol, val = gl.swarm_optimize(env)
    # out = gl.apply_rotatron_solution(sol, env, protein.copy())

    print("Optimization complete, writing output PDB...")

    out.to_pdb("/Users/noahhk/Downloads/fold_mev_nag_model_0_protein_glycosylated.pdb")
    print("Glycosylation complete.")

    out.save("/Users/noahhk/Downloads/fold_mev_nag_model_0_protein_glycosylated.pkl")


def test_glycanshield_symelas_protein():
    out = gl.Protein.load(
        "/Users/noahhk/Downloads/fold_mev_nag_model_0_protein_glycosylated.pkl"
    )

    quickshield = gl.quickshield(
        out,
        repeats=3,
        angle_step=120,
        save_conformations_to="/Users/noahhk/Downloads/fold_mev_nag_model_0_protein_glycosylated_conformations",
    )
    print("Glycan shield simulation complete.")
    quickshield.df.to_csv(
        "/Users/noahhk/Downloads/fold_mev_nag_model_0_protein_glycosylated_exposure.csv"
    )
