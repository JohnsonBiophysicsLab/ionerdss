# Geometry Conventions

This note defines the signed angles `phi` and `omega` used in ioNERDSS bond geometry.

## Setup

For a binding pair, define

- `v1 = bind_site1 - com1`
- `v2 = bind_site2 - com2`
- `sigma1 = bind_site1 - bind_site2`
- `sigma2 = -sigma1`
- `n1`, `n2` as the local normal vectors for molecule 1 and molecule 2

The vectors `v1` and `v2` point from each molecule center to its interface site. The vector `sigma1` points from interface 2 to interface 1.

## What is `phi`?

`phi1` and `phi2` describe how each molecule's local normal is rotated about its own interface vector.

For molecule 1, define

```math
t^{(\mathrm{plane})}_1 = \frac{v_1 \times \sigma_1}{\lVert v_1 \times \sigma_1 \rVert},
\qquad
t^{(\mathrm{norm})}_1 = \frac{v_1 \times n_1}{\lVert v_1 \times n_1 \rVert}.
```

Then the unsigned angle is

```math
\phi_1 = \cos^{-1}\!\left(t^{(\mathrm{plane})}_1 \cdot t^{(\mathrm{norm})}_1\right).
```

Analogously, for molecule 2,

```math
t^{(\mathrm{plane})}_2 = \frac{v_2 \times \sigma_2}{\lVert v_2 \times \sigma_2 \rVert},
\qquad
t^{(\mathrm{norm})}_2 = \frac{v_2 \times n_2}{\lVert v_2 \times n_2 \rVert},
```

```math
\phi_2 = \cos^{-1}\!\left(t^{(\mathrm{plane})}_2 \cdot t^{(\mathrm{norm})}_2\right).
```

Interpretation: `phi` is the rotation from the binding plane `(v, sigma)` to the local orientation plane `(v, n)`, measured around the axis `v`.

## What is `omega`?

`omega` is the torsion angle between the two interface vectors `v1` and `v2`, measured around the inter-site axis `sigma1`.

Define

```math
a_1 = \frac{\sigma_1 \times v_1}{\lVert \sigma_1 \times v_1 \rVert},
\qquad
a_2 = \frac{\sigma_1 \times v_2}{\lVert \sigma_1 \times v_2 \rVert}.
```

Then the unsigned torsion is

```math
\omega = \cos^{-1}\!\left(a_1 \cdot a_2\right).
```

Interpretation: `omega` measures the twist needed to rotate the projected `v1` direction into the projected `v2` direction about the axis `sigma1`.

## Torsion Direction Convention

The sign convention is right-handed and axis-based.

For `phi1`, project `n1` and `sigma1` onto the plane perpendicular to `v1`:

```math
n_{1,\perp} = n_1 - (n_1 \cdot \hat v_1)\hat v_1,
\qquad
\sigma_{1,\perp} = \sigma_1 - (\sigma_1 \cdot \hat v_1)\hat v_1.
```

Then define the direction test

```math
d_{\phi_1} = \frac{\sigma_{1,\perp} \times n_{1,\perp}}
{\lVert \sigma_{1,\perp} \times n_{1,\perp} \rVert}.
```

If `d_{phi1}` is parallel to `v1`, ioNERDSS stores `phi1` as negative. If it is anti-parallel to `v1`, `phi1` stays positive. The same construction is used for `phi2` with `(v2, sigma2, n2)`.

For `omega`, project `v1` and `v2` onto the plane perpendicular to `sigma1`:

```math
v_{1,\perp} = v_1 - (v_1 \cdot \hat \sigma_1)\hat \sigma_1,
\qquad
v_{2,\perp} = v_2 - (v_2 \cdot \hat \sigma_1)\hat \sigma_1.
```

Then define

```math
d_\omega = \frac{v_{1,\perp} \times v_{2,\perp}}
{\lVert v_{1,\perp} \times v_{2,\perp} \rVert}.
```

If `d_omega` is parallel to `sigma1`, ioNERDSS stores `omega` as negative. If it is anti-parallel to `sigma1`, `omega` stays positive.

So, in this codebase:

- positive `omega` means the right-hand normal of `(v1_perp, v2_perp)` points opposite to `sigma1`
- negative `omega` means that normal points along `sigma1`

This is the implemented convention, even if another text might choose the opposite sign.

## Range and Degenerate Cases

- `phi1`, `phi2`, and `omega` are stored in radians.
- The underlying unsigned angle is in `[0, \pi]`; the sign rule extends this to a signed convention.
- If a cross product used in the construction vanishes, the torsion is geometrically degenerate.
- In some linear-molecule cases, `phi` may be undefined and represented as `NaN`.

## Relation to the Implementation

These conventions are implemented in:

- `ionerdss/utils/bond_geometry.py`
- `ionerdss/model/pdb/nerdss_exporter.py`

Those functions should be treated as the source of truth for exported ioNERDSS geometry.
