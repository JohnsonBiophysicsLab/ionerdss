# NERDSS Pipeline Troubleshooting Guide

This guide addresses common warnings and errors encountered during the ionerdss pipeline execution.

## Warnings

##### `Representative instance {instance_id} missing interfaces: {missing_interfaces}. Falling back to mixed-frame instances (RISKY).`

**Cause:**
The pipeline attempts to align all instances of a molecule to a single "representative" reference frame to define relative interface positions. However, the chosen representative instance lacks some interfaces present in other instances (e.g., due to disorder or missing density in the PDB).

**Implication:**
The coordinates for the missing interfaces will be derived from *other* instances that do have them via kabsch transformation. This can lead to inconsistencies if the internal geometry of the molecule varies significantly between instances (e.g., flexible domains). It's "RISKY" because the relative positions might not be perfectly rigorous.

**Solution/Action:**
1.  **Ignore if Stable:** If the molecule is rigid, this warning is often benign.
2.  **Input Quality:** Ensure your input PDB/CIF is complete and high quality.

---

##### `small sigma values : {sigma_value}; consider increasing binding radius threshold`

**Cause:**
The calculated standard deviation (`sigma`) of the bond length for an interface is very small (< 0.5 nm). This usually happens when:
1.  There are very few instances of this bond (statistics are poor).
2.  The binding sites are extremely constrained/clashed in the PDB.
3.  The `binding_radius_threshold` used for detection is too small, cutting off the distribution.

**Implication:**
In ioNERDSS, a very small `sigma` makes the bond geometry extremely sensitive. e.g. a small shift in the binding site position in the PDB can lead to a large shift in all the related angles and torsion angles.

**Solution/Action:**
1.  **Increase Threshold:** Increase `interface_detect_distance_cutoff` in your hyperparameters (e.g., from 1.0 to 1.5 or 2.0 nm). (and `interface_detect_n_residue_cutoff` accordingly)
2. **Check Clashes:** Visualize the structure to ensure the interface isn't physically impossible (clashed). If so, turn on the steric clash mode. `steric_clash_mode="auto"`

---

##### `No matched interfaces to define rotation for {chain_id}`

**Cause:**
The code is trying to align a chain to the reference frame but cannot find enough shared interfaces to calculate a uniquely defined rotation matrix. This might happen for a chain that only has 1 inerface, or the interfaces are colinear.

**Implication:**
The code is going to pick an arbitrary rotation matrix satisfying the existing shared interfaces. This might lead to a wrong relative position of the chain in the complex.

**Solution/Action:**
1.  **Check Connectivity:** Ensure this chain is actually part of the complex you intend to simulate and only has 1 interfaces or interfaces colinear with the center of mass.
2.  **Check Missing Interfaces:** Sometimes this may occur due to missing interfaces. Increasing `interface_detect_distance_cutoff` might help find the missing interfaces.

---

##### `ODE pipeline skipped (continuing with normal workflow): Assembly has {num_molecules} molecules, exceeding max_complex_size ({max_complex_size}). Skipping ODE generation...`

**Cause:**
The automated ODE generation (which builds a reaction network) converts the system graph into a reaction system. This process scales exponentially with complex size. To prevent hanging, it has a safety limit (`max_complex_size`, default 12).

**Implication:**
The simulation *will* run (NERDSS files are generated), but you won't get the automated ODE system and solution for this specific system.

**Solution/Action:**
1.  **Increase Limit (Risky):** If you really need it, increase `max_complex_size` in the `generate_ode_model_from_system` call, but be aware it might take a long time or crash.
2.  **Use Simulation:** Rely on the NERDSS simulation (which handles large complexes natively) instead of the ODE approximation.

## Errors
