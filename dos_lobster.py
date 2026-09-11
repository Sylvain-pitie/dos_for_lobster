#!/usr/bin/env python3
"""
dos_plot.py -- Trace la DOS totale + projetée (élément/orbitale) à partir de
DOSCAR.lobster, avec la même mise en page que cobi_plot.py (énergie en
ordonnée) et la charte de couleurs des DOS HSE06 du SI.

Charte par défaut :
    TDOS  -> noir
    H  s  -> rouge
    Na s  -> jaune
    Na p  -> orange
    Ae s  -> bleu ciel      (Ae = Be, Mg, Ca, Sr, Ba)
    Ae p  -> bleu foncé
    Ae d  -> violet         (seulement si présent dans la base LOBSTER)

Fichiers requis dans le répertoire courant : DOSCAR.lobster et POSCAR.
"""

import argparse
import os
import re
import sys

import numpy as np
from matplotlib import pyplot as plt

from pymatgen.core import Structure
from pymatgen.io.lobster import Doscar
from pymatgen.electronic_structure.dos import Dos
from pymatgen.electronic_structure.core import Spin, OrbitalType

# --------------------------------------------------------------------------
# Charte graphique (éditable)
# --------------------------------------------------------------------------
AE_ELEMENTS = {"Be", "Mg", "Ca", "Sr", "Ba"}

COLORS = {
    "Total": "#000000",   # noir
    "H_s":   "#E31A1C",   # rouge
    "Na_s":  "#FFC20A",   # jaune
    "Na_p":  "#FF7F00",   # orange
    "Ae_s":  "#56B4E9",   # bleu ciel
    "Ae_p":  "#00308F",   # bleu foncé
    "Ae_d":  "#7B3294",   # violet
}

FALLBACK_COLORS = ["#4DAF4A", "#984EA3", "#A65628", "#999999", "#F781BF"]

ORBITAL_TYPES = {"s": OrbitalType.s, "p": OrbitalType.p,
                 "d": OrbitalType.d, "f": OrbitalType.f}


def color_for(el, orb, counter=[0]):
    """Couleur d'une projection (élément, orbitale)."""
    key = f"{el}_{orb}"
    if key in COLORS:
        return COLORS[key]
    if el in AE_ELEMENTS and f"Ae_{orb}" in COLORS:
        return COLORS[f"Ae_{orb}"]
    c = FALLBACK_COLORS[counter[0] % len(FALLBACK_COLORS)]
    counter[0] += 1
    return c


# --------------------------------------------------------------------------
# Utilitaires
# --------------------------------------------------------------------------
def element_orbital_types(cdos):
    """
    Parcourt les pDOS LOBSTER et renvoie
    {element: {'s': {n, ...}, 'p': {n, ...}}} où n est le nombre quantique
    principal effectivement présent dans la base LOBSTER.
    """
    out = {}
    for site, orb_dict in cdos.pdos.items():
        try:
            el = site.specie.symbol
        except AttributeError:
            el = site.species_string
        for orb in orb_dict:
            name = str(orb).strip()
            m = re.match(r"^(\d+)\s*[-_]?\s*([spdfSPDF])", name)
            if m:
                out.setdefault(el, {}).setdefault(m.group(2).lower(), set()).add(int(m.group(1)))
                continue
            # Repli : nom sans nombre quantique principal (ex. 'p_x', 'px')
            m = re.search(r"([spdf])", name.lower())
            if m:
                out.setdefault(el, {}).setdefault(m.group(1), set())
    return out


def orbital_label(el, orb, ns):
    """'Na', 'p', {2} -> 'Na 2p'."""
    if ns:
        return f"{el} {'/'.join(str(n) for n in sorted(ns))}{orb}"
    return f"{el} {orb}"


def ordered_elements(available):
    """H et Na en tête, puis les alcalino-terreux, puis le reste."""
    order = []
    for el in ("H", "Na"):
        if el in available:
            order.append(el)
    for el in sorted(available):
        if el not in order and el in AE_ELEMENTS:
            order.append(el)
    for el in sorted(available):
        if el not in order:
            order.append(el)
    return order


def parse_selection(spec, available):
    """
    'auto' -> toutes les orbitales présentes dans la base.
    'H:s;Na:s,p;Be:s,p' -> sélection explicite.
    Renvoie [(element, orbitale), ...]
    """
    selection = []
    if spec.strip().lower() == "auto":
        for el in ordered_elements(available):
            for orb in ("s", "p", "d", "f"):
                if orb in available[el]:
                    selection.append((el, orb))
        return selection

    for block in spec.split(";"):
        block = block.strip()
        if not block:
            continue
        if ":" not in block:
            raise ValueError(f"Bloc mal formé dans --orbitals : '{block}' "
                             "(attendu 'El:s,p')")
        el, orbs = block.split(":", 1)
        el = el.strip()
        if el not in available:
            print(f"  Attention : élément {el} absent de DOSCAR.lobster, ignoré")
            continue
        for orb in orbs.split(","):
            orb = orb.strip().lower()
            if orb not in ORBITAL_TYPES:
                print(f"  Attention : orbitale '{orb}' inconnue, ignorée")
                continue
            if orb not in available[el]:
                print(f"  Attention : {el} {orb} absent de la base LOBSTER, ignoré")
                continue
            selection.append((el, orb))
    return selection


def band_edges(energies, dens_total, e0=0.0, tol=1e-3):
    """VBM / CBM détectés de part et d'autre de e0 avec un seuil sur la DOS."""
    dens_total = np.asarray(dens_total, dtype=float)
    vbm = cbm = None
    idx = np.where((energies <= e0) & (dens_total > tol))[0]
    if len(idx):
        vbm = float(energies[idx].max())
    idx = np.where((energies >= e0) & (dens_total > tol))[0]
    if len(idx):
        cbm = float(energies[idx].min())
    return vbm, cbm


def orbital_letter(name):
    """'2p_x' -> 'p'."""
    m = re.match(r"^(\d+)\s*[-_]?\s*([spdfSPDF])", str(name).strip())
    if m:
        return m.group(2).lower()
    m = re.search(r"([spdf])", str(name).lower())
    return m.group(1) if m else None


def group_site_dos(cdos, sites_ordered, indices, orbs=None, average=False):
    """
    Somme les pDOS d'une liste de sites (indices 1-based, numérotation de
    --list-orbitals), éventuellement restreinte à certains types d'orbitales.
    """
    total = {}
    n = 0
    for i in indices:
        if i < 1 or i > len(sites_ordered):
            raise ValueError(f"Indice de site hors limites : {i} "
                             f"(1..{len(sites_ordered)})")
        site = sites_ordered[i - 1]
        n += 1
        for orb, dens in cdos.pdos[site].items():
            if orbs and orbital_letter(orb) not in orbs:
                continue
            for spin, arr in dens.items():
                arr = np.asarray(arr, dtype=float)
                total[spin] = total[spin] + arr if spin in total else arr.copy()
    if average and n:
        total = {s: d / n for s, d in total.items()}
    return total, n


def classify_sites(structure, target, neighbour, cutoff):
    """
    Regroupe les sites de l'élément `target` selon leur nombre de voisins
    `neighbour` dans un rayon `cutoff`. Renvoie {n_voisins: [indices 1-based]}.
    """
    classes = {}
    for i, site in enumerate(structure, start=1):
        try:
            el = site.specie.symbol
        except AttributeError:
            el = site.species_string
        if el != target:
            continue
        count = sum(1 for nb in structure.get_neighbors(site, cutoff)
                    if nb.specie.symbol == neighbour)
        classes.setdefault(count, []).append(i)
    return classes


def save_command(output_filename, command_line):
    cmd_filename = output_filename + ".command"
    with open(cmd_filename, "w") as f:
        f.write(command_line + "\n")
    print(f"Commande sauvegardée dans {cmd_filename}")


# --------------------------------------------------------------------------
# Programme principal
# --------------------------------------------------------------------------
def plot_dos(args, command_line):
    print(f"Lecture de {args.doscar} et {args.poscar}...")
    doscar = Doscar(doscar=args.doscar, structure_file=args.poscar)
    cdos = doscar.completedos

    energies = np.array(cdos.energies, dtype=float)
    efermi = float(cdos.efermi)

    print(f"  E_F lue dans l'en-tête : {efermi:.4f} eV")
    print(f"  Grille d'énergie brute : [{energies.min():.3f}, {energies.max():.3f}] eV")

    # --- Référencement à E_F -------------------------------------------------
    if args.eref == "efermi":
        shift = efermi
    elif args.eref == "none":
        shift = 0.0
    else:  # auto
        shift = efermi if abs(efermi) > 0.1 else 0.0
    if shift != 0.0:
        print(f"  -> énergies décalées de -{shift:.4f} eV (E_F = 0)")
    else:
        print("  -> pas de décalage (DOSCAR déjà référencé à E_F = 0)")
    energies = energies - shift

    # --- Repositionnement du zéro -------------------------------------------
    dens_total = sum(np.asarray(d, dtype=float) for d in cdos.densities.values())
    vbm, cbm = band_edges(energies, dens_total, e0=0.0, tol=args.dos_tol)
    if vbm is not None and cbm is not None:
        print(f"  VBM = {vbm:.4f} eV, CBM = {cbm:.4f} eV, "
              f"gap = {cbm - vbm:.4f} eV (seuil DOS = {args.dos_tol})")

    extra = 0.0
    if args.zero == "vbm":
        if vbm is None:
            sys.exit("VBM non détecté : ajuste --dos-tol ou utilise --zero efermi")
        extra = vbm
    elif args.zero == "cbm":
        if cbm is None:
            sys.exit("CBM non détecté : ajuste --dos-tol ou utilise --zero efermi")
        extra = cbm
    elif args.zero == "gapcenter":
        if vbm is None or cbm is None:
            sys.exit("Bords de bande non détectés : utilise --zero efermi")
        extra = 0.5 * (vbm + cbm)

    extra += args.eshift
    if extra != 0.0:
        print(f"  -> zéro repositionné sur '{args.zero}' : décalage "
              f"supplémentaire de {extra:.4f} eV")
        energies = energies - extra
    print(f"  Décalage total depuis la grille brute : {shift + extra:.4f} eV")
    print(f"  (pour aligner cobi_plot.py : --eshift {shift + extra:.4f})")

    # --- Sélection des projections ------------------------------------------
    available = element_orbital_types(cdos)

    if args.list_orbitals:
        print("\nNoms d'orbitales bruts lus par pymatgen, site par site :")
        for i, (site, orb_dict) in enumerate(cdos.pdos.items(), start=1):
            try:
                el = site.specie.symbol
            except AttributeError:
                el = site.species_string
            print(f"  site {i:>3}  {el:<3} : {sorted(str(o) for o in orb_dict)}")
        print("\nBase reconstruite :",
              {el: {o: sorted(ns) for o, ns in d.items()} for el, d in available.items()})
        return

    print("  Base LOBSTER détectée :",
          {el: {o: sorted(ns) for o, ns in d.items()} for el, d in available.items()})

    selection = parse_selection(args.orbitals, available)
    if not selection:
        print("Aucune projection sélectionnée.")

    # --- Construction des courbes -------------------------------------------
    curves = []  # (label, couleur, {Spin: densités}, is_total)

    if not args.no_total:
        dens = (cdos.get_smeared_densities(args.sigma) if args.sigma
                else cdos.densities)
        curves.append(("TDOS", COLORS["Total"], dens, True))

    for el, orb in selection:
        spd = cdos.get_element_spd_dos(el)
        otype = ORBITAL_TYPES[orb]
        if otype not in spd:
            print(f"  Attention : {el} {orb} non renvoyé par get_element_spd_dos, ignoré")
            continue
        dos_obj = spd[otype]
        dens = (dos_obj.get_smeared_densities(args.sigma) if args.sigma
                else dos_obj.densities)
        label = orbital_label(el, orb, available[el].get(orb, set()))
        curves.append((label, color_for(el, orb), dens, False))

    # --- Courbes résolues par site ------------------------------------------
    sites_ordered = list(cdos.pdos.keys())
    structure = getattr(doscar, "structure", None)
    if structure is None:
        structure = getattr(cdos, "structure", None)
    if structure is None:
        structure = Structure.from_file(args.poscar)

    # Cohérence entre l'ordre des sites des pDOS et celui de la structure
    if len(structure) != len(sites_ordered):
        sys.exit(f"Incohérence : {len(structure)} sites dans {args.poscar} "
                 f"mais {len(sites_ordered)} dans les pDOS")

    site_specs = []
    if args.classify:
        bits = args.classify.split(":")
        if len(bits) < 2:
            sys.exit("--classify attend 'VOISIN:CUTOFF[:CIBLE]', ex. 'Be:1.8:H'")
        neighbour, cutoff = bits[0].strip(), float(bits[1])
        target = bits[2].strip() if len(bits) > 2 else "H"
        classes = classify_sites(structure, target, neighbour, cutoff)

        class_names = {}
        if args.classify_names:
            for item in args.classify_names.split(";"):
                if "=" not in item:
                    continue
                k, v = item.split("=", 1)
                class_names[int(k.strip())] = v.strip()

        class_colors = {}
        if args.classify_colors:
            for item in args.classify_colors.split(";"):
                if "=" not in item:
                    continue
                k, v = item.split("=", 1)
                class_colors[int(k.strip())] = v.strip()

        print(f"\n  Classification des {target} par nombre de {neighbour} "
              f"dans {cutoff} Å :")
        for count in sorted(classes):
            idx = classes[count]
            name = class_names.get(count, f"{target} ({count} {neighbour})")
            print(f"    {count} {neighbour} -> {len(idx)} sites : {idx}  [{name}]")
            site_specs.append((name, idx, None, class_colors.get(count)))
        if not classes:
            print("    aucun site trouvé")

    for block in (args.sites.split(";") if args.sites else []):
        block = block.strip()
        if not block:
            continue
        bits = block.split(":")
        if len(bits) < 2:
            sys.exit(f"Bloc --sites mal formé : '{block}' "
                     "(attendu 'Nom:indices[:orbitales][:couleur]')")
        name = bits[0].strip()
        indices = []
        for part in bits[1].split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-")
                indices.extend(range(int(a), int(b) + 1))
            elif part:
                indices.append(int(part))
        orbs = ([o.strip().lower() for o in bits[2].split(",") if o.strip()]
                if len(bits) > 2 and bits[2].strip() else None)
        col = bits[3].strip() if len(bits) > 3 and bits[3].strip() else None
        site_specs.append((name, indices, orbs, col))

    for k, (name, indices, orbs, col) in enumerate(site_specs):
        dens, nsites = group_site_dos(cdos, sites_ordered, indices, orbs,
                                      average=args.site_average)
        if not dens:
            print(f"  Attention : groupe '{name}' vide, ignoré")
            continue
        if args.sigma:
            dens = Dos(cdos.efermi, energies, dens).get_smeared_densities(args.sigma)
        color = col or FALLBACK_COLORS[k % len(FALLBACK_COLORS)]
        label = f"{name}" + (" /at." if args.site_average else "")
        print(f"  Groupe '{name}' : {nsites} sites"
              + (f", orbitales {orbs}" if orbs else ""))
        curves.append((label, color, dens, False))

    if args.names:
        names = [n.strip() for n in args.names.split(";")]
        if len(names) != len(curves):
            raise ValueError(f"--names contient {len(names)} noms pour "
                             f"{len(curves)} courbes")
        curves = [(n, c, d, t) for n, (_, c, d, t) in zip(names, curves)]

    spin_polarized = any(Spin.down in d for _, _, d, _ in curves)
    print(f"  Calcul {'spin-polarisé' if spin_polarized else 'non spin-polarisé'}, "
          f"{len(curves)} courbe(s)")

    spin_mode = args.spin_mode
    if spin_mode == "auto":
        spin_mode = "mirror" if spin_polarized else "sum"

    # --- Tracé ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=tuple(args.figsize))

    mask = np.ones_like(energies, dtype=bool)
    if args.ymin is not None:
        mask &= energies >= args.ymin
    if args.ymax is not None:
        mask &= energies <= args.ymax

    all_x = []
    for label, color, dens, is_total in curves:
        if spin_mode == "sum":
            channels = [(sum(dens[s] for s in dens), "-", label)]
        elif spin_mode == "mirror":
            channels = [(dens[Spin.up], "-", label)]
            if Spin.down in dens:
                channels.append((-np.asarray(dens[Spin.down]), "-", None))
        else:  # updown
            channels = [(dens[Spin.up], "-", label)]
            if Spin.down in dens:
                channels.append((dens[Spin.down], "--", None))

        for y, ls, lab in channels:
            y = np.asarray(y, dtype=float)
            ax.plot(y, energies, color=color, linestyle=ls,
                    linewidth=args.linewidth, label=lab, zorder=3 if is_total else 4)
            if is_total and args.fill:
                ax.fill_betweenx(energies, 0, y, color=color, alpha=0.10, zorder=2)
            if np.any(mask):
                all_x.extend(y[mask])

    if str(args.zero_line).strip().lower() != "none":
        ax.axhline(y=float(args.zero_line), color="black", linestyle="dashed",
                   linewidth=1.5, zorder=1)
    if spin_mode == "mirror" and str(args.xmin).strip().lower() == "auto":
        ax.axvline(x=0, color="black", linestyle="dashed", linewidth=1.5, zorder=1)

    ax.legend(frameon=False, fontsize=args.fontsize, loc=args.legend_loc)
    ax.tick_params(labelsize=args.fontsize)
    ax.set_ylabel(args.ylabel, fontsize=args.fontsize)
    ax.set_xlabel(args.xlabel, fontsize=args.fontsize)

    if args.ymin is not None or args.ymax is not None:
        ax.set_ylim(args.ymin, args.ymax)

    auto_xmin = str(args.xmin).strip().lower() == "auto"
    if auto_xmin:
        xmin = min(min(all_x), 0.0) if all_x else 0.0
    else:
        xmin = float(args.xmin)

    xmax = args.dosmax if args.dosmax is not None else (max(all_x) if all_x else None)

    if xmax is not None:
        span = (xmax - xmin) or 1.0
        left = xmin - 0.05 * span if auto_xmin else xmin
        right = xmax if args.dosmax is not None else xmax + 0.05 * span
        ax.set_xlim(left, right)

    fig.subplots_adjust(top=0.98, bottom=0.10, left=0.14, right=0.96)

    out = args.output + ".png"
    fig.savefig(out, format="png", dpi=args.dpi)
    if args.pdf:
        fig.savefig(args.output + ".pdf", format="pdf")
        print(f"Figure vectorielle sauvegardée sous {args.output}.pdf")
    print(f"Figure sauvegardée sous {out}")
    if args.show:
        plt.show()

    save_command(args.output, command_line)


def main():
    p = argparse.ArgumentParser(
        description="Trace la DOS totale et projetée à partir de DOSCAR.lobster.")
    p.add_argument("--doscar", default="DOSCAR.lobster")
    p.add_argument("--poscar", default="POSCAR")
    p.add_argument("--orbitals", default="auto",
                   help="'auto' ou sélection explicite, ex. 'H:s;Na:s,p;Be:s,p'")
    p.add_argument("--names", type=str, default=None,
                   help="Noms des courbes séparés par ';' (TDOS incluse si tracée)")
    p.add_argument("--sites", type=str, default=None,
                   help="Groupes de sites, séparés par ';'. Format "
                        "'Nom:indices[:orbitales][:couleur]', indices 1-based "
                        "(numérotation de --list-orbitals), plages avec '-'. "
                        "Ex. 'H bridge:25-30:s:#E31A1C;H term:31-36:s:#FB9A99'")
    p.add_argument("--classify", type=str, default=None,
                   help="Classe automatiquement des sites par nombre de voisins : "
                        "'VOISIN:CUTOFF[:CIBLE]', ex. 'Be:1.8:H' regroupe les H "
                        "selon leur nombre de Be à moins de 1.8 Å")
    p.add_argument("--classify-names", type=str, default=None,
                   help="Renomme les classes de --classify. Format "
                        "'n=Nom' séparés par ';'. Spécifique à la structure. "
                        "Ex. '1=H (term);2=H (bridge)'")
    p.add_argument("--classify-colors", type=str, default=None,
                   help="Couleurs des classes de --classify. Format 'n=couleur' "
                        "séparés par ';'. Ex. '1=#E31A1C;2=#FB9A99'")
    p.add_argument("--site-average", action="store_true",
                   help="Normalise les groupes de sites par atome")
    p.add_argument("--list-orbitals", action="store_true",
                   help="Affiche les noms d'orbitales lus dans DOSCAR.lobster et quitte")
    p.add_argument("--no-total", action="store_true", help="Ne pas tracer la TDOS")
    p.add_argument("--fill", action="store_true",
                   help="Remplissage léger sous la TDOS")
    p.add_argument("--sigma", type=float, default=None,
                   help="Élargissement gaussien en eV (ex. 0.05)")
    p.add_argument("--eref", choices=["auto", "efermi", "none"], default="auto",
                   help="Traitement de l'E_F lue dans l'en-tête (défaut: auto)")
    p.add_argument("--zero", choices=["efermi", "vbm", "cbm", "gapcenter"],
                   default="efermi",
                   help="Où placer le zéro d'énergie (défaut: efermi)")
    p.add_argument("--eshift", type=float, default=0.0,
                   help="Décalage manuel supplémentaire en eV (E -> E - eshift)")
    p.add_argument("--zero-line", default="0",
                   help="Ordonnée de la ligne pointillée de référence, "
                        "ou 'none' pour la masquer (défaut: 0)")
    p.add_argument("--dos-tol", type=float, default=1e-3,
                   help="Seuil de DOS pour la détection VBM/CBM (défaut: 1e-3)")
    p.add_argument("--spin-mode", choices=["auto", "sum", "mirror", "updown"],
                   default="auto")
    p.add_argument("--ymin", type=float, default=None)
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--xmin", default="0",
                   help="Borne gauche de l'axe DOS (défaut: 0, collée à l'axe ; "
                        "'auto' pour laisser apparaître les valeurs négatives)")
    p.add_argument("--dosmax", type=float, default=None,
                   help="Limite de l'axe DOS (défaut: automatique sur la fenêtre)")
    p.add_argument("--fontsize", type=int, default=20)
    p.add_argument("--linewidth", type=float, default=3.0)
    p.add_argument("--figsize", type=float, nargs=2, default=[6, 8])
    p.add_argument("--legend-loc", default="best")
    p.add_argument("--xlabel", default="DOS")
    p.add_argument("--ylabel", default="$E - E_\\mathrm{F}$ (eV)")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--pdf", action="store_true", help="Sortie PDF vectorielle en plus")
    p.add_argument("--show", action="store_true")
    p.add_argument("--output", default="testdos")
    args = p.parse_args()

    for f in (args.doscar, args.poscar):
        if not os.path.isfile(f):
            sys.exit(f"Fichier introuvable : {f}")

    plot_dos(args, " ".join(sys.argv))


if __name__ == "__main__":
    main()
