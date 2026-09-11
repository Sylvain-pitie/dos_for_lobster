# dos_plot.py

Plot total and projected densities of states from LOBSTER output, with energy on
the vertical axis and a colour scheme designed for publication figures.

Beyond the usual element- and orbital-projected DOS, the script can resolve the
DOS by **crystallographic site**, including automatic grouping of atoms by their
local coordination — for example separating bridging from terminal hydrogen.

---

## Contents

- [Requirements](#requirements)
- [Input files](#input-files)
- [Quick start](#quick-start)
- [Selecting projections](#selecting-projections)
- [Site-resolved DOS](#site-resolved-dos)
- [Energy reference](#energy-reference)
- [Spin-polarised calculations](#spin-polarised-calculations)
- [Colours](#colours)
- [Negative DOS values](#negative-dos-values)
- [Troubleshooting](#troubleshooting)
- [Option reference](#option-reference)

---

## Requirements

```
python >= 3.8
pymatgen
numpy
matplotlib
```

```bash
pip install pymatgen numpy matplotlib
chmod +x dos_plot.py
```

Tested with LOBSTER 5.1.1 and VASP 6.x.

---

## Input files

Run from the directory holding your LOBSTER output.

| File | Purpose |
|------|---------|
| `DOSCAR.lobster` | Projected DOS on the local basis |
| `POSCAR` | Structure: site ordering and neighbour analysis |

The `POSCAR` must be the one used for the LOBSTER run. A different species
ordering silently shifts the site-to-element mapping.

Every run writes `<output>.png` and `<output>.command`, the latter containing the
full command line for reproducibility.

---

## Quick start

Always start by inspecting what the file actually contains:

```bash
dos_plot.py --list-orbitals
```

```
  site   1  Na  : ['2p_x', '2p_y', '2p_z', '3s']
  site  13  Be  : ['2s', '2p_x', '2p_y', '2p_z']
  site  25  H   : ['1s']

Reconstructed basis: {'Na': {'s': [3], 'p': [2]}, 'Be': {'s': [2], 'p': [2]}, 'H': {'s': [1]}}
```

This tells you which projections exist and gives the site numbering used
elsewhere. Then plot:

```bash
dos_plot.py --orbitals 'H:s;Na:s,p;Be:s,p' \
            --ymin -6 --ymax 8 --fontsize 14 --linewidth 2 \
            --output dos_NaBeH3
```

---

## Selecting projections

`--orbitals auto` (the default) plots every orbital type present in the basis.
To choose explicitly:

```bash
--orbitals 'H:s;Na:s,p;Be:s,p'
```

Legend labels carry the real principal quantum numbers read from the file, so
you get `Na 3s`, `Na 2p`, `Be 2s` rather than a generic `Na s`. Override every
label with `--names`, one per curve, in plotting order, including the total DOS
if it is plotted:

```bash
--names 'TDOS;H 1s;Na 3s;Na 2p;Be 2s;Be 2p'
```

`--no-total` omits the total DOS; `--fill` adds a light fill beneath it;
`--sigma 0.05` applies Gaussian broadening in eV.

Note that the total DOS written by LOBSTER is the sum of the projections onto
the local basis, not the VASP total DOS. The two differ by the charge spilling.

---

## Site-resolved DOS

Two ways to separate chemically distinct atoms of the same element.

### Automatic, by coordination

```bash
--classify 'Be:1.8:H'
```

reads as *group the H sites by how many Be lie within 1.8 Å*. The syntax is
`NEIGHBOUR:CUTOFF[:TARGET]`, `TARGET` defaulting to `H`. Periodic boundary
conditions are respected, so an atom bridging a cell boundary is handled
correctly.

The classes are reported with their site indices:

```
  H classification by number of Be within 1.8 A:
    1 Be -> 6 sites : [31, 32, 33, 34, 35, 36]  [H (term)]
    2 Be -> 6 sites : [25, 26, 27, 28, 29, 30]  [H (bridge)]
```

Choose the cutoff from your tabulated bond lengths: too short and everything
becomes terminal, too long and everything becomes bridging.

Class names and colours are **structure-specific**, so they are not hard-coded.
Map them yourself, keyed by the neighbour count:

```bash
--classify-names '1=H (term);2=H (bridge)' \
--classify-colors '1=#E31A1C;2=#FB9A99'
```

A class you do not name keeps its automatic label — a useful warning that the
cutoff produced an unexpected category, such as a third class of H.

### Manual, by site index

```bash
--sites 'H (bridge):25-30:s:#FB9A99;H (term):31-36:s:#E31A1C'
```

Format: `Name:indices[:orbitals][:colour]`, groups separated by `;`. Indices are
1-based, follow the `--list-orbitals` numbering, and accept ranges with `-`.
Orbitals and colour are optional. `--sites` and `--orbitals` combine freely, so
you can keep Na and Be as element projections while detailing hydrogen.

### Normalisation

`--site-average` divides each site group by its number of atoms. Use it whenever
groups differ in size, otherwise you compare intensities that scale with atom
count rather than with local electronic structure.

---

## Energy reference

Three independent controls, applied in order.

**1. `--eref {auto,efermi,none}`** handles the Fermi energy in the
`DOSCAR.lobster` header. `auto` shifts only if it differs from zero by more than
0.1 eV. LOBSTER usually writes energies already referenced to E_F, so `auto`
normally does nothing. The script prints what it did and the raw energy range —
check this on the first run of a new system.

**2. `--zero {efermi,vbm,cbm,gapcenter}`** moves the zero onto a band edge. For
an insulator, E_F sits at an arbitrary position inside the gap, so `--zero vbm`
is often the meaningful choice. Edges are located by thresholding the total DOS:

```
  VBM = -0.0501 eV, CBM = 5.2566 eV, gap = 5.3066 eV (DOS threshold = 0.001)
```

The detected gap depends on `--dos-tol`. Before quoting a value, confirm it is
stable across a decade of thresholds (`1e-4`, `1e-3`, `1e-2`); a numerical tail
counted as a real state will narrow the gap artificially.

**3. `--eshift X`** applies an extra manual shift, `E -> E - X`.

`--zero-line` places the dashed horizontal reference line anywhere on the final
scale, or hides it with `none`.

The script prints the cumulative shift, ready to paste into a companion COBI
plot so both panels share one energy scale:

```
  total shift from raw grid: -0.0501 eV
  (to align cobi_plot.py: --eshift -0.0501)
```

---

## Spin-polarised calculations

`--spin-mode auto` mirrors the spin-down channel below zero when the calculation
is spin-polarised, and sums the channels otherwise.

| Mode | Behaviour |
|------|-----------|
| `sum` | Adds up and down, everything positive |
| `updown` | Both positive, spin-down dashed |
| `mirror` | Spin-down drawn negative |

`--xmin` defaults to `0`, which clips everything below zero. With `mirror` you
must pass `--xmin auto`, or the spin-down channel falls outside the frame.

For a closed-shell compound the calculation should not be spin-polarised at all.
If the console reports one, `ISPIN = 2` was set and both channels are identical:
prefer `updown` or rerun with `ISPIN = 1`, since `sum` would double the values.

---

## Colours

Defaults live in the `COLORS` dictionary at the top of the file:

| Curve | Colour |
|-------|--------|
| Total DOS | black |
| H s | red |
| Na s | yellow |
| Na p | orange |
| Ae s (Be, Mg, Ca, Sr, Ba) | light blue |
| Ae p | dark blue |
| Ae d | purple |

`Ae_s`, `Ae_p` and `Ae_d` are generic keys covering every alkaline-earth element,
so one command works across a whole series without editing. Site groups take
colours from `FALLBACK_COLORS` unless set through `--classify-colors` or the
fourth field of `--sites`.

When plotting two subsets of one element, two shades of that element's colour
keep the figure readable and consistent with the rest of a paper. Avoid plotting
the global element projection alongside its own subsets: it is exactly their sum
and adds nothing but a third line to distinguish.

---

## Negative DOS values

Small negative values in the *projected* DOS are expected. LOBSTER projects onto
a non-orthogonal local basis and partitions the overlap Mulliken-style, which
can yield negative local contributions near band edges or in antibonding
regions. As long as the amplitude stays at the level of a few percent of the
maximum, this is the usual artefact of the partitioning.

The *total* DOS is never negative. If it is, you are looking at a mirrored
spin-down channel — see the spin section above.

`--xmin 0` (the default) clips negative values from view. This is cosmetic
only: the curves are unchanged, and `--xmin auto` restores the full range if you
want to judge the amplitude first.

---

## Troubleshooting

**An orbital is missing from the projection.** Run `--list-orbitals`; it prints
the orbital names per site exactly as they appear in the file. If an expected
orbital is absent, the calculation is at fault, not the plotting.

**`basisfunctions` in `lobsterin` seem ignored.** If `projectionData.lobster`
exists and `loadProjectionFromFile` is set, LOBSTER reloads the stored
projection and never applies your current basis lines. Delete or rename the
file, or comment out the keyword, then rerun. Keep only `saveProjectionToFile`
until the basis is settled. Verify what was used:

```bash
grep -A20 -i "basis functions" lobsterout
```

**An unoccupied orbital cannot be added.** Tabulated basis sets provide
functions only for orbitals occupied in the free atom, so Ba 5d is unavailable
(Ba is [Xe]6s²). This mainly affects the conduction band, but matters under
pressure where s→d transfer populates those states.

**`AttributeError: 'Doscar' object has no attribute 'structure'`.** Older
pymatgen versions do not expose it; the script falls back to reading `POSCAR`
directly. Make sure the file is present and matches the run.

**Site indices look wrong.** They follow the POSCAR ordering. The script checks
that the site count matches the pDOS and stops if not.

---

## Option reference

| Option | Default | Description |
|--------|---------|-------------|
| `--doscar` | `DOSCAR.lobster` | Input DOS file |
| `--poscar` | `POSCAR` | Structure file |
| `--orbitals` | `auto` | `auto`, or `'H:s;Na:s,p'` |
| `--names` | – | Rename all curves, `;`-separated, in plot order |
| `--sites` | – | `Name:indices[:orbitals][:colour]`, `;`-separated |
| `--classify` | – | `NEIGHBOUR:CUTOFF[:TARGET]`, e.g. `Be:1.8:H` |
| `--classify-names` | – | `n=Name`, `;`-separated |
| `--classify-colors` | – | `n=colour`, `;`-separated |
| `--site-average` | off | Normalise site groups per atom |
| `--list-orbitals` | off | Print orbital names per site and exit |
| `--no-total` | off | Omit the total DOS |
| `--fill` | off | Light fill under the total DOS |
| `--sigma` | – | Gaussian broadening in eV |
| `--eref` | `auto` | `auto`, `efermi`, `none` |
| `--zero` | `efermi` | `efermi`, `vbm`, `cbm`, `gapcenter` |
| `--eshift` | `0.0` | Extra manual shift, `E -> E - X` |
| `--zero-line` | `0` | Reference line position, or `none` |
| `--dos-tol` | `1e-3` | DOS threshold for band-edge detection |
| `--spin-mode` | `auto` | `auto`, `sum`, `mirror`, `updown` |
| `--ymin`, `--ymax` | – | Energy window in eV |
| `--xmin` | `0` | Left DOS bound, or `auto` |
| `--dosmax` | – | Right DOS bound |
| `--fontsize` | `20` | Axes and legend font size |
| `--linewidth` | `3.0` | Curve width |
| `--figsize` | `6 8` | Figure size in inches |
| `--legend-loc` | `best` | Matplotlib legend location |
| `--xlabel` | `DOS` | x-axis label |
| `--ylabel` | `$E - E_\mathrm{F}$ (eV)` | y-axis label |
| `--dpi` | `300` | PNG resolution |
| `--pdf` | off | Also write a vector PDF |
| `--show` | off | Open an interactive window |
| `--output` | `testdos` | Output basename |

---

## Citing

- R. Nelson, C. Ertural, J. George, V. L. Deringer, G. Hautier, R. Dronskowski,
  *J. Comput. Chem.* **2020**, 41, 1931.
- S. P. Ong et al., *Comput. Mater. Sci.* **2013**, 68, 314.
