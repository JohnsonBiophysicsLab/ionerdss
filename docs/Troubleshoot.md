# NERDSS Pipeline Troubleshooting Guide

This guide addresses common warnings and errors encountered during the ionerdss pipeline execution.

## Warnings

##### `Representative instance {instance_id} missing interfaces: {missing_interfaces}. Falling back to mixed-frame instances (RISKY).`

**Cause:**
The pipeline attempts to align all instances of a molecule to a single "representative" reference frame to define relative interface positions. However, the chosen representative instance lacks some interfaces present in other instances (e.g., due to disorder or missing density in the PDB).

**Implication:**
The coordinates for the missing interfaces will be derived from *other* instances that do have them, rotated into the representative's frame by a Kabsch fit of the interfaces the two instances share (or, when they share fewer than two, by the rotation between the chains' Cα-aligned frames). This can lead to inconsistencies if the internal geometry of the molecule varies significantly between instances (e.g., flexible domains). It's "RISKY" because the relative positions might not be perfectly rigorous.

**Solution/Action:**
1.  **Ignore if Stable:** If the molecule is rigid, this warning is often benign.
2.  **Input Quality:** Ensure your input PDB/CIF is complete and high quality.

---

##### `small sigma values : {sigma_value}; consider increasing binding radius threshold`

**Cause:**
The calculated bond length (`sigma`, the distance between the two binding sites, averaged over all bound pairs of that reaction) is very small (< 0.5 nm). This usually happens when:
1.  The binding sites are extremely constrained/clashed in the PDB.
2.  The `interface_detect_distance_cutoff` used for detection is small, so each binding site sits at the centroid of only the few residues closest to the partner.

**Implication:**
In ioNERDSS, a very small `sigma` makes the bond geometry extremely sensitive. e.g. a small shift in the binding site position in the PDB can lead to a large shift in all the related angles and torsion angles.

**Solution/Action:**
1.  **Increase Threshold:** Increase `interface_detect_distance_cutoff` in your hyperparameters (e.g., from the default 0.9 to 1.5 or 2.0 nm). (and `interface_detect_n_residue_cutoff` accordingly)
2. **Check Clashes:** Visualize the structure to ensure the interface isn't physically impossible (clashed). If so, turn on the steric clash mode. `steric_clash_mode="auto"`

---

##### `No matched interfaces to define rotation for {molecule_name}`

**Cause:**
The code is trying to align a chain to the reference frame of its molecule's representative instance. It fits the interfaces the two share, and falls back to the rotation between the chains' Cα-aligned frames when they share fewer than two. This warning means the chain shares no interface type with the representative and that fallback is unavailable (e.g. the chain's Cα alignment failed).

**Implication:**
No rotation is applied: the chain's binding-site normal vectors fall back to the representative's unrotated ones, and any interface borrowed from it for the representative (see the mixed-frame warning above) keeps its unrotated local coordinates. This might lead to a wrong relative position of the chain in the complex.

**Solution/Action:**
1.  **Check Connectivity:** Ensure this chain is actually part of the complex you intend to simulate.
2.  **Check Missing Interfaces:** Sometimes this may occur due to missing interfaces. Increasing `interface_detect_distance_cutoff` might help find the missing interfaces.

---

##### `ODE pipeline skipped (continuing with normal workflow): Assembly has {num_molecules} molecules, exceeding max_complex_size ({max_complex_size}). Skipping ODE generation...`

**Cause:**
The automated ODE generation (which builds a reaction network) converts the system graph into a reaction system. This process scales exponentially with complex size. To prevent hanging, it has a safety limit (the `max_complex_size_ode` hyperparameter, default 12).

**Implication:**
The simulation *will* run (NERDSS files are generated), but you won't get the automated ODE system and solution for this specific system.

**Solution/Action:**
1.  **Increase Limit (Risky):** If you really need it, increase `max_complex_size_ode` in your hyperparameters (or `max_complex_size` if you call `generate_ode_model_from_system` yourself), but be aware it might take a long time or crash.
2.  **Use Simulation:** Rely on the NERDSS simulation (which handles large complexes natively) instead of the ODE approximation.

## Errors
