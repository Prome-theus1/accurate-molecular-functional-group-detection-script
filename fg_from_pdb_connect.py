#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import csv
from collections import defaultdict, Counter
from dataclasses import dataclass
from pathlib import Path

from typing import Dict, List, Tuple, Optional, Set, Any

# A fixed, stable set of FG keys to always include in CSV outputs (even if count = 0)
BASE_FG_KEYS: List[str] = [
    # Rings (5-8), mutually exclusive categories
    # 5/6-membered rings use 3 classes: heteroring (any hetero atom), aromatic ring (no hetero + aromatic), aliphatic ring (no hetero + non-aromatic)
    "heteroring5", "aromatic_ring5", "aliphatic_ring5",
    "heteroring6", "aromatic_ring6", "aliphatic_ring6",
    "ring7",
    "ring8",

    # O-related
    "phenol", "alcohol", "ether", "epoxide",

    # N-related
    "nitrile", "nitro", "amide", "amine",

    # S-related
    "thioether", "sulfone_or_sulfonyl", "sulfoxide_or_sulfenyl_oxide", "thioester",

    # P-related
    "phosphate", "phosphorus_containing",

    # Halogens
    "halide",

    # Carbonyl-centered
    "carboxylic_acid", "ester", "acyl_chloride",
    "ketone", "aldehyde",
]


# ----------------------------
# Basic PDB parsing
# ----------------------------

@dataclass
class Atom:
    serial: int
    name: str
    element: str
    chain: str
    resseq: int
    resname: str
    x: float
    y: float
    z: float


def _safe_int(s: str, default: int = 0) -> int:
    try:
        return int(s)
    except Exception:
        return default


def _safe_float(s: str, default: float = float("nan")) -> float:
    try:
        return float(s)
    except Exception:
        return default


def guess_element(atom_name: str, element_field: str) -> str:
    """
    Prefer PDB element column (77-78). If empty, guess from atom name.
    """
    el = (element_field or "").strip()
    if el:
        return el.capitalize()

    # Guess from atom name: e.g. " O1 " -> O, "CL1" -> Cl
    name = (atom_name or "").strip()
    if not name:
        return "X"
    # Two-letter elements
    two = name[:2].strip().capitalize()
    if two in {"Cl", "Br", "Si", "Na", "Ca", "Li", "Al", "Mg", "Zn", "Fe", "Cu", "Mn", "Co", "Ni", "Se"}:
        return two
    # One-letter
    one = name[0].upper()
    return one


def parse_pdb_atoms_and_conect(pdb_path: Path, ignore_h: bool = True):
    """
    Returns:
      atoms: dict serial -> Atom (H optionally removed)
      edges_raw: list of (i, j) from CONECT (may include duplicates)
      serials_in_file_order: list[int] (all ATOM/HETATM serials in file order, including H if ignore_h=False)
    """
    atoms: Dict[int, Atom] = {}
    edges_raw: List[Tuple[int, int]] = []
    serials_in_file_order: List[int] = []

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                serial = _safe_int(line[6:11].strip())
                name = line[12:16]
                resname = line[17:20].strip()
                chain = (line[21:22] or "").strip() or "_"
                resseq = _safe_int(line[22:26].strip(), 0)
                x = _safe_float(line[30:38].strip())
                y = _safe_float(line[38:46].strip())
                z = _safe_float(line[46:54].strip())
                element_field = line[76:78] if len(line) >= 78 else ""
                element = guess_element(name, element_field)

                serials_in_file_order.append(serial)
                if ignore_h and element.upper() == "H":
                    continue

                atoms[serial] = Atom(
                    serial=serial,
                    name=name.strip(),
                    element=element.capitalize(),
                    chain=chain,
                    resseq=resseq,
                    resname=resname,
                    x=x, y=y, z=z
                )

            elif line.startswith("CONECT"):
                # Format: CONECT a b c d e ...
                parts = line.split()
                if len(parts) >= 3:
                    a = _safe_int(parts[1], -1)
                    for p in parts[2:]:
                        b = _safe_int(p, -1)
                        if a > 0 and b > 0:
                            edges_raw.append((a, b))

    return atoms, edges_raw, serials_in_file_order


# ----------------------------
# Optional OpenBabel bond order + aromaticity + rings
# ----------------------------

def get_openbabel_bond_info(
    pdb_path: Path,
    serials_in_file_order: List[int],
    atoms_all: Dict[int, Atom],
):
    """
    Use OpenBabel to perceive bond orders and aromatic atoms from a PDB.

    Returns
    -------
    bond_order : dict[(min_serial, max_serial)] -> 1/2/3
    aromatic_atoms : set[serial]
    []  # rings_info is always empty; ring detection is not handled here.
    """
    try:
        from openbabel import openbabel as ob
    except Exception:
        return {}, set(), []

    conv = ob.OBConversion()
    if not conv.SetInFormat("pdb"):
        return {}, set(), []

    mol = ob.OBMol()
    ok = conv.ReadFile(mol, str(pdb_path))
    if not ok:
        return {}, set(), []

    # perceive
    mol.PerceiveBondOrders()
    mol.FindSSSR()

    # Aromaticity perception differs across OpenBabel Python bindings.
    try:
        mol.PerceiveAromaticity()
    except AttributeError:
        try:
            ob.PerceiveAromaticity(mol)
        except Exception:
            pass

    # --- robust index->serial mapping ---
    # Sometimes OpenBabel may drop explicit hydrogens; try both with-H and no-H mapping.
    serials_with_h = list(serials_in_file_order)
    serials_no_h = [s for s in serials_in_file_order if (s in atoms_all and atoms_all[s].element.upper() != "H")]

    chosen_serials: Optional[List[int]] = None
    if mol.NumAtoms() == len(serials_with_h):
        chosen_serials = serials_with_h
    elif mol.NumAtoms() == len(serials_no_h):
        chosen_serials = serials_no_h
    else:
        # cannot safely map
        return {}, set(), []

    idx_to_serial = {idx: chosen_serials[idx - 1] for idx in range(1, mol.NumAtoms() + 1)}

    bond_order: Dict[Tuple[int, int], int] = {}
    aromatic_atoms: Set[int] = set()

    # aromatic atoms
    for atom in ob.OBMolAtomIter(mol):
        if atom.IsAromatic():
            aromatic_atoms.add(idx_to_serial[atom.GetIdx()])

    # bonds
    for bond in ob.OBMolBondIter(mol):
        i = idx_to_serial[bond.GetBeginAtomIdx()]
        j = idx_to_serial[bond.GetEndAtomIdx()]
        a, b = (i, j) if i < j else (j, i)
        bo = int(bond.GetBondOrder())
        # OB sometimes returns 5 for aromatic; clamp to 1 for our FG rules
        if bo not in (1, 2, 3):
            bo = 1
        bond_order[(a, b)] = bo

    # rings_info is always empty; ring detection is not handled here.
    return bond_order, aromatic_atoms, []


def dist(a: Atom, b: Atom) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


# ----------------------------
# Bond order estimation (fallback heuristic)
# ----------------------------

def estimate_bond_order(el1: str, el2: str, d: float, bo_override: Optional[int] = None) -> int:
    """
    Very rough bond order estimation from bond length (Å).
    Works best for C=O, N=O, C#N, etc.
    Returns 1,2,3.
    """
    if bo_override in (1, 2, 3):
        return int(bo_override)

    a, b = el1.upper(), el2.upper()
    pair = tuple(sorted([a, b]))

    if pair == ("C", "O"):
        return 2 if d < 1.30 else 1

    if pair == ("C", "N"):
        if d < 1.20:
            return 3
        if d < 1.32:
            return 2
        return 1

    if pair == ("N", "O"):
        return 2 if d < 1.25 else 1

    if pair == ("O", "S"):
        return 2 if d < 1.52 else 1

    return 1


# ----------------------------
# Functional group rules (graph + local environment)
# ----------------------------

def residue_key(atom: Atom):
    return (atom.chain, atom.resseq, atom.resname)


def build_graph(atoms: Dict[int, Atom], edges_raw: List[Tuple[int, int]]):
    """
    Build adjacency with edge multiplicity from raw CONECT.
    Returns:
      adj: dict i -> set(j)
      mult: dict (min(i,j),max(i,j)) -> count
    """
    adj: Dict[int, Set[int]] = defaultdict(set)
    mult: Counter = Counter()

    for i, j in edges_raw:
        if i not in atoms or j not in atoms:
            continue
        if i == j:
            continue
        a, b = (i, j) if i < j else (j, i)
        mult[(a, b)] += 1
        adj[i].add(j)
        adj[j].add(i)

    return adj, mult


# ----------------------------
# Ring detection from CONECT graph (no OpenBabel)
# ----------------------------

def _heavy_adj(atoms: Dict[int, Atom], adj: Dict[int, Set[int]]) -> Dict[int, Set[int]]:
    """Adjacency restricted to heavy atoms (exclude H)."""
    heavy = {i for i, a in atoms.items() if a.element.upper() != "H"}
    hadj: Dict[int, Set[int]] = defaultdict(set)
    for i in heavy:
        for j in adj.get(i, set()):
            if j in heavy:
                hadj[i].add(j)
    return hadj


def _cycle_edges_from_path(path: List[int]) -> frozenset:
    """Given a path that starts and ends at same node, build an undirected edge-set key."""
    edges = []
    for u, v in zip(path[:-1], path[1:]):
        a, b = (u, v) if u < v else (v, u)
        edges.append((a, b))
    return frozenset(edges)


def find_simple_cycles_upto_n(
    atoms: Dict[int, Atom],
    adj: Dict[int, Set[int]],
    max_len: int = 8,
) -> List[Set[int]]:
    """Find unique simple cycles (undirected) up to length max_len using DFS.

    This is a lightweight molecular-graph ring finder based on CONECT.
    It avoids duplicates by:
      - Only allowing cycles where the start node is the smallest serial in that cycle.
      - Canonicalizing cycles by their undirected edge-set.

    Returns a list of cycles as sets of atom serials (heavy atoms only).
    """
    hadj = _heavy_adj(atoms, adj)
    nodes = sorted(hadj.keys())

    seen_edge_sets: Set[frozenset] = set()
    cycles: List[Set[int]] = []

    for start in nodes:
        # DFS stack: (current, parent, path)
        stack: List[Tuple[int, int, List[int]]] = [(start, -1, [start])]
        while stack:
            cur, parent, path = stack.pop()
            if len(path) > max_len:
                continue

            for nb in hadj.get(cur, set()):
                if nb == parent:
                    continue

                # Ensure start is the smallest node in the cycle to reduce duplicates
                if nb < start:
                    continue

                if nb == start:
                    # close a cycle
                    if 3 <= len(path) <= max_len:
                        cyc_path = path + [start]
                        edge_key = _cycle_edges_from_path(cyc_path)
                        if edge_key not in seen_edge_sets:
                            seen_edge_sets.add(edge_key)
                            node_set = set(path)
                            cycles.append(node_set)
                    continue

                if nb in path:
                    continue

                stack.append((nb, cur, path + [nb]))

    return cycles


def ring_aromatic_heuristic(
    atoms: Dict[int, Atom],
    hadj: Dict[int, Set[int]],
    ring_nodes: Set[int],
) -> bool:
    """Very rough aromaticity heuristic without bond orders.

    It is intentionally conservative and only used to populate aromatic_* ring keys.
    If you don't trust it, just rely on non-aromatic ring keys.
    """
    size = len(ring_nodes)
    if size not in (5, 6, 7, 8):
        return False

    # Basic element + degree sanity
    deg3 = 0
    for i in ring_nodes:
        el = atoms[i].element.upper()
        if el == "H":
            return False
        if el not in ("C", "N", "O", "S", "P"):
            return False
        d = len(hadj.get(i, set()))
        if d < 2 or d > 3:
            return False
        if d == 3:
            deg3 += 1

    # Fused/aryl rings tend to have more degree-3 atoms
    if size == 6:
        return deg3 >= 3
    if size == 5:
        return deg3 >= 2
    if size == 7:
        return deg3 >= 4
    if size == 8:
        return deg3 >= 5
    return False


def detect_rings_from_graph(
    atoms: Dict[int, Atom],
    adj: Dict[int, Set[int]],
    aromatic_atoms: Optional[Set[int]] = None,
) -> List[Dict[str, Any]]:
    """Detect rings (sizes 5-8) from CONECT-only graph.

    Returns rings_info list of dicts with:
      {'size': int, 'aromatic': bool, 'hetero': bool, 'serials': set[int]}

    Each ring contributes to exactly one category later (ring*/heteroring*/aromatic_*).
    """
    hadj = _heavy_adj(atoms, adj)
    cycles = find_simple_cycles_upto_n(atoms, adj, max_len=8)
    aromatic_atoms = aromatic_atoms or set()

    rings_info: List[Dict[str, Any]] = []
    for cyc in cycles:
        size = len(cyc)
        if size not in (5, 6, 7, 8):
            continue
        hetero = any(atoms[i].element.upper() not in ("C", "H") for i in cyc)
        if aromatic_atoms:
            # Prefer OpenBabel aromatic-atom marks: if a strict majority of atoms in the ring are aromatic, label the ring aromatic.
            arom_n = sum(1 for i in cyc if i in aromatic_atoms)
            aromatic = (arom_n * 2 > size)
        else:
            aromatic = ring_aromatic_heuristic(atoms, hadj, cyc)
        rings_info.append({
            "size": size,
            "aromatic": aromatic,
            "hetero": hetero,
            "serials": set(cyc),
        })

    return rings_info


def classify_residue(
    atoms: Dict[int, Atom],
    adj: Dict[int, Set[int]],
    bond_order_map: Optional[Dict[Tuple[int, int], int]] = None,
    aromatic_atoms: Optional[Set[int]] = None,
    rings_info: Optional[List[Dict[str, Any]]] = None,
):
    """
    For a residue, scan non-C atoms and infer FG types.
    Returns Counter of FG occurrences.
    """
    fg: Counter = Counter()

    bond_order_map = bond_order_map or {}
    aromatic_atoms = aromatic_atoms or set()
    rings_info = rings_info or []

    # Precompute ring membership for this residue for robust aromatic-carbon decisions
    ring_atoms_all: Set[int] = set()
    aromatic_ring_atoms_all: Set[int] = set()
    for r in rings_info:
        s = set(r.get("serials", set()) or set())
        if not s:
            continue
        ring_atoms_all |= s
        if bool(r.get("aromatic", False)):
            aromatic_ring_atoms_all |= s

    def get_bo(i: int, j: int) -> Optional[int]:
        a, b = (i, j) if i < j else (j, i)
        return bond_order_map.get((a, b))

    # Pre-index atoms by residue
    atoms_by_res: Dict[Tuple[str, int, str], List[int]] = defaultdict(list)
    for aid, a in atoms.items():
        atoms_by_res[residue_key(a)].append(aid)

    # Helper: find carbonyl carbons in residue (C with O double)
    def is_carbonyl_carbon(c_id: int) -> bool:
        c = atoms[c_id]
        if c.element.upper() != "C":
            return False
        for nb in adj.get(c_id, set()):
            b = atoms[nb]
            if b.element.upper() == "O":
                bo = get_bo(c_id, nb)
                bo = estimate_bond_order("C", "O", dist(c, b), bo_override=bo)
                if bo == 2:
                    return True
        return False

    # Helper: list double-bonded oxygens on a carbon
    def carbonyl_oxygens(c_id: int) -> List[int]:
        res: List[int] = []
        c = atoms[c_id]
        for nb in adj.get(c_id, set()):
            b = atoms[nb]
            if b.element.upper() == "O":
                bo = get_bo(c_id, nb)
                bo = estimate_bond_order("C", "O", dist(c, b), bo_override=bo)
                if bo == 2:
                    res.append(nb)
        return res

    # Aromatic carbon: require ring membership to avoid false positives (e.g., ethanol)
    def is_probably_aromatic_carbon(c_id: int) -> bool:
        c = atoms[c_id]
        if c.element.upper() != "C":
            return False

        # Must be part of a detected ring to be considered aromatic
        if c_id not in ring_atoms_all:
            return False

        # If OpenBabel marked it aromatic, accept (but only within a ring)
        if c_id in aromatic_atoms:
            return True

        # Otherwise, use our ring-based heuristic aromatic labeling
        if c_id in aromatic_ring_atoms_all:
            return True

        # Last resort: heavy-atom degree within the residue graph
        heavy_deg = sum(1 for nb in adj.get(c_id, set()) if atoms[nb].element.upper() != "H")
        return heavy_deg >= 3

    def ring_key(size: int, aromatic: bool, hetero: bool) -> Optional[str]:
        """Make ring classification mutually exclusive.

        For 5/6-membered rings, use 3 classes:
          1) aromatic_ring{size}: aromatic regardless of hetero atoms
          2) heteroring{size}: non-aromatic + hetero
          3) aliphatic_ring{size}: non-aromatic + no hetero

        For 7/8-membered rings, keep the legacy collapsed keys: ring7 / ring8.
        """
        if size not in (5, 6, 7, 8):
            return None

        # For 7/8 membered rings, collapse everything into ring7/ring8
        if size in (7, 8):
            return f"ring{size}"

        # For 5/6 membered rings, enforce 3-class scheme with priority:
        #   1) aromatic_ring{size}: aromatic regardless of hetero atoms
        #   2) heteroring{size}: non-aromatic + hetero
        #   3) aliphatic_ring{size}: non-aromatic + no hetero
        if aromatic:
            return f"aromatic_ring{size}"
        if hetero:
            return f"heteroring{size}"
        return f"aliphatic_ring{size}"

    for _resk, atom_ids in atoms_by_res.items():
        phenol_carbons: Set[int] = set()  # aromatic carbon(s) bonded to phenolic OH in this residue
        # Scan hetero atoms first (H is NOT a center atom; only used as neighbor)
        for aid in atom_ids:
            a = atoms[aid]
            el = a.element.upper()
            if el in ("C", "H"):
                continue

            nbs = list(adj.get(aid, set()))
            nb_elems = [atoms[x].element.upper() for x in nbs]
            nb_count = Counter(nb_elems)

            # ------ Oxygen rules ------
            if el == "O":
                c_neighbors = [x for x in nbs if atoms[x].element.upper() == "C"]

                # O single to C: alcohol/phenol/ether-like oxygen
                if len(c_neighbors) == 1:
                    c_id = c_neighbors[0]
                    bo0 = get_bo(c_id, aid)
                    bo = estimate_bond_order("C", "O", dist(atoms[c_id], a), bo_override=bo0)
                    if bo == 2:
                        # carbonyl oxygen handled in carbonyl-centered classification
                        pass
                    else:
                        # If this oxygen is single-bonded to a carbonyl carbon (e.g., carboxylic acid / ester),
                        # skip O-centered classification to avoid double counting; carbonyl-centered block handles it.
                        if is_carbonyl_carbon(c_id):
                            continue
                        has_h = (nb_count.get("H", 0) >= 1)
                        if is_probably_aromatic_carbon(c_id):
                            if has_h:
                                fg["phenol"] += 1
                                phenol_carbons.add(c_id)
                            else:
                                # aryl_ether category removed by request
                                pass
                        else:
                            if has_h:
                                fg["alcohol"] += 1
                            else:
                                fg["ether"] += 1

                # O connected to two carbons: epoxide vs ether (mutually exclusive)
                if nb_count["C"] == 2 and len(c_neighbors) == 2:
                    # Avoid counting carbonyl-adjacent O (ester/carbonyl systems)
                    if any(is_carbonyl_carbon(c) for c in c_neighbors):
                        pass
                    else:
                        c1, c2 = c_neighbors
                        is_epoxide = (c2 in adj.get(c1, set()))
                        if is_epoxide:
                            fg["epoxide"] += 1
                        else:
                            # Do not count cyclic (ring) oxygens as ether; ring category already captures them
                            if aid not in ring_atoms_all:
                                fg["ether"] += 1

                # NOTE: 删除 sulfur_oxygen_bond / N_O_bond 的计数（按你要求）

            # ------ Nitrogen rules ------
            elif el == "N":
                is_nitrile_n = False
                is_amide_n = False
                is_nitro_n = False

                # nitrile: N triple to C
                for nb in nbs:
                    b = atoms[nb]
                    if b.element.upper() == "C":
                        bo0 = get_bo(nb, aid)
                        bo = estimate_bond_order("C", "N", dist(b, a), bo_override=bo0)
                        if bo == 3:
                            is_nitrile_n = True
                            fg["nitrile"] += 1

                # nitro: N connected to >=2 O
                if nb_count["O"] >= 2:
                    is_nitro_n = True
                    fg["nitro"] += 1

                # amide: N single-bonded to carbonyl carbon
                for nb in nbs:
                    if atoms[nb].element.upper() == "C" and is_carbonyl_carbon(nb):
                        is_amide_n = True
                        fg["amide"] += 1

                # amine: N connected to >=1 C AND not amide AND not nitrile AND not nitro AND not in any ring
                if (nb_count["C"] >= 1) and (not is_amide_n) and (not is_nitrile_n) and (not is_nitro_n) and (aid not in ring_atoms_all):
                    fg["amine"] += 1

            # ------ Sulfur rules ------
            elif el == "S":
                # Mutually exclusive sulfur categories (priority: sulfone > sulfoxide > thioether)
                # thiol category removed by request
                o_n = int(nb_count.get("O", 0))
                c_n = int(nb_count.get("C", 0))

                if o_n >= 2:
                    fg["sulfone_or_sulfonyl"] += 1
                elif o_n == 1:
                    fg["sulfoxide_or_sulfenyl_oxide"] += 1
                else:
                    # thioether: S bonded to >=2 carbons, but exclude thioesters (S bonded to a carbonyl carbon)
                    if c_n >= 2 and (aid not in ring_atoms_all):
                        c_neighbors = [x for x in nbs if atoms[x].element.upper() == "C"]
                        if any(is_carbonyl_carbon(c) for c in c_neighbors):
                            pass
                        else:
                            fg["thioether"] += 1
                # NOTE: 删除 sulfur_mono_carbon 的计数（按你要求）

            # ------ Phosphorus rules ------
            elif el == "P":
                if nb_count["O"] >= 3:
                    fg["phosphate"] += 1
                else:
                    fg["phosphorus_containing"] += 1

            # ------ Halogens ------
            elif el in ("F", "CL", "BR", "I"):
                # Count any halogen as `halide` regardless of bonding,
                # except acyl chlorides (Cl bonded to a carbonyl carbon), which are handled as `acyl_chloride`.
                if el == "CL":
                    c_neighbors = [x for x in nbs if atoms[x].element.upper() == "C"]
                    if any(is_carbonyl_carbon(c) for c in c_neighbors):
                        pass
                    else:
                        fg["halide"] += 1
                else:
                    fg["halide"] += 1

            else:
                fg[f"element_{el}"] += 1

        # Ring features (from CONECT graph): count rings fully contained in this residue
        # If this residue contains phenol, do NOT additionally count the benzene ring as aromatic_ring6.
        atom_set = set(atom_ids)
        for r in rings_info:
            serials = r.get("serials", set())
            if not serials or (not serials.issubset(atom_set)):
                continue
            size = int(r.get("size", 0))
            aromatic = bool(r.get("aromatic", False))
            hetero = bool(r.get("hetero", False))

            # Skip aromatic benzene-ring counting for phenol
            if size == 6 and aromatic and phenol_carbons:
                if any(c in serials for c in phenol_carbons):
                    continue

            k = ring_key(size, aromatic, hetero)
            if k is not None:
                fg[k] += 1

        # ---- Carbonyl-centered classification ----
        for cid in atom_ids:
            c = atoms[cid]
            if c.element.upper() != "C":
                continue
            o_dbl = carbonyl_oxygens(cid)
            if not o_dbl:
                continue

            neighbors = [nb for nb in adj.get(cid, set()) if nb not in o_dbl]
            nb_e = [atoms[nb].element.upper() for nb in neighbors]
            nb_count = Counter(nb_e)

            # acid / ester
            if nb_count["O"] >= 1:
                o_single = [nb for nb in neighbors if atoms[nb].element.upper() == "O"]
                if o_single:
                    for o in o_single:
                        o_nbs = [x for x in adj.get(o, set()) if x != cid]
                        if any(atoms[x].element.upper() == "C" for x in o_nbs):
                            fg["ester"] += 1
                        else:
                            fg["carboxylic_acid"] += 1

            # (amide counting via carbonyl C removed to avoid double-counting)

            if nb_count["CL"] >= 1:
                fg["acyl_chloride"] += 1

            if nb_count["S"] >= 1:
                fg["thioester"] += 1

            # Refined ketone/aldehyde: only if NOT acyl derivative (no single-bond O/N/S/Cl)
            is_acyl_derivative = (nb_count["O"] >= 1) or (nb_count["N"] >= 1) or (nb_count["S"] >= 1) or (nb_count["CL"] >= 1)

            c_neighbors = int(nb_count["C"])
            h_neighbors = int(nb_count.get("H", 0))

            if not is_acyl_derivative:
                if c_neighbors >= 2 and h_neighbors == 0:
                    fg["ketone"] += 1
                elif c_neighbors == 1 and h_neighbors == 1:
                    fg["aldehyde"] += 1
                else:
                    # carbonyl_other removed by request
                    pass

    return fg


# ----------------------------
# Batch helper
# ----------------------------

def infer_fgs_for_pdb(pdb_path: Path, use_openbabel: bool):
    """Infer functional groups for a single PDB and return a Counter of FG counts (summed over residues)."""
    atoms, edges_raw, serials_in_file_order = parse_pdb_atoms_and_conect(pdb_path, ignore_h=False)

    bond_order_map, aromatic_atoms, _ = ({}, set(), [])
    if use_openbabel:
        bond_order_map, aromatic_atoms, _ = get_openbabel_bond_info(
            pdb_path, serials_in_file_order, atoms
        )

    adj, _mult = build_graph(atoms, edges_raw)
    rings_info = detect_rings_from_graph(atoms, adj, aromatic_atoms=aromatic_atoms)

    res_atoms = defaultdict(list)
    for aid, a in atoms.items():
        res_atoms[residue_key(a)].append(aid)

    fg_total = Counter()
    for (_chain, _resseq, _resname), atom_ids in res_atoms.items():
        sub_atoms = {aid: atoms[aid] for aid in atom_ids}
        sub_adj = defaultdict(set)
        for aid in atom_ids:
            for nb in adj.get(aid, set()):
                if nb in sub_atoms:
                    sub_adj[aid].add(nb)

        fg = classify_residue(
            sub_atoms,
            sub_adj,
            bond_order_map=bond_order_map,
            aromatic_atoms=aromatic_atoms,
            rings_info=rings_info,
        )
        fg_total.update(fg)

    ob_ok = bool(use_openbabel and bond_order_map and len(bond_order_map) > 0)
    return fg_total, ob_ok


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Infer functional groups from PDB using CONECT by local bond environment (H kept as auxiliary neighbors)."
    )
    ap.add_argument("--pdb", default=None, help="Input PDB file (with CONECT).")
    ap.add_argument("--pdb-dir", default=None, help="Input directory containing PDB files to batch process.")
    ap.add_argument("--recursive", action="store_true", help="Recursively search for *.pdb under --pdb-dir.")

    ap.add_argument("--out-annotated", default="annotated.csv", help="(Single file) Per-residue annotation CSV.")
    ap.add_argument("--out-summary", default="summary.csv", help="(Single file) Summary stats output CSV.")

    ap.add_argument("--out-table", default="fg_table.csv",
                    help="(Batch mode) Output CSV table: one row per PDB, columns are FG counts.")
    ap.add_argument("--out-total", default="fg_total_summary.csv",
                    help="(Batch mode) Output CSV summary: total FG counts across all PDBs.")

    ap.add_argument(
        "--use-openbabel",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Use OpenBabel to perceive bond orders/aromatic atoms. Ring detection uses CONECT graph; ring aromaticity uses OB aromatic atoms when available.",
    )

    args = ap.parse_args()

    if (args.pdb is None) and (args.pdb_dir is None):
        raise ValueError("Please provide either --pdb (single file) or --pdb-dir (batch mode).")
    if (args.pdb is not None) and (args.pdb_dir is not None):
        raise ValueError("Please provide only one of --pdb or --pdb-dir, not both.")

    # ---------------- Batch mode ----------------
    if args.pdb_dir is not None:
        pdb_dir = Path(args.pdb_dir)
        if not pdb_dir.exists():
            raise ValueError(f"--pdb-dir not found: {pdb_dir}")

        pdb_files = sorted(pdb_dir.rglob("*.pdb") if args.recursive else pdb_dir.glob("*.pdb"))
        if not pdb_files:
            raise ValueError(f"No *.pdb files found under: {pdb_dir}")

        per_file = []  # list of (relpath, Counter)
        total = Counter()
        all_keys = set()
        ob_on_count = 0

        for p in pdb_files:
            fg_counter, ob_ok = infer_fgs_for_pdb(p, args.use_openbabel)
            rel = str(p.relative_to(pdb_dir))
            per_file.append((rel, fg_counter))
            total.update(fg_counter)
            all_keys.update(fg_counter.keys())
            if args.use_openbabel and ob_ok:
                ob_on_count += 1

        # Always include base FG keys, then extras
        extra_keys = sorted([k for k in all_keys if k not in BASE_FG_KEYS])
        fg_keys = list(BASE_FG_KEYS) + extra_keys
        header = ["pdb"] + fg_keys
        with open(args.out_table, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            for rel, c in per_file:
                row = {"pdb": rel}
                for k in fg_keys:
                    row[k] = int(c.get(k, 0))
                w.writerow(row)

        # write total summary, all fg_keys in order, zeros explicit
        with open(args.out_total, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["functional_group", "count"])
            w.writeheader()
            for k in fg_keys:
                w.writerow({"functional_group": k, "count": int(total.get(k, 0))})

        print(f"[OK] batch pdb files: {len(pdb_files)}")
        print(f"[OK] openbabel: {'ON' if args.use_openbabel else 'OFF'} (ok for {ob_on_count}/{len(pdb_files)})")
        print(f"[OK] wrote: {args.out_table}")
        print(f"[OK] wrote: {args.out_total}")
        return

    # ---------------- Single-file mode ----------------
    pdb_path = Path(args.pdb)
    atoms, edges_raw, serials_in_file_order = parse_pdb_atoms_and_conect(pdb_path, ignore_h=False)

    bond_order_map, aromatic_atoms, _ = ({}, set(), [])
    if args.use_openbabel:
        bond_order_map, aromatic_atoms, _ = get_openbabel_bond_info(
            pdb_path, serials_in_file_order, atoms
        )

    adj, _mult = build_graph(atoms, edges_raw)
    rings_info = detect_rings_from_graph(atoms, adj, aromatic_atoms=aromatic_atoms)

    # group atoms by residue
    res_atoms = defaultdict(list)
    for aid, a in atoms.items():
        res_atoms[residue_key(a)].append(aid)

    all_summary = Counter()
    annotated_rows = []

    for (chain, resseq, resname), atom_ids in sorted(res_atoms.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        sub_atoms = {aid: atoms[aid] for aid in atom_ids}
        sub_adj = defaultdict(set)
        for aid in atom_ids:
            for nb in adj.get(aid, set()):
                if nb in sub_atoms:
                    sub_adj[aid].add(nb)

        fg = classify_residue(
            sub_atoms,
            sub_adj,
            bond_order_map=bond_order_map,
            aromatic_atoms=aromatic_atoms,
            rings_info=rings_info,
        )

        if not fg:
            continue

        all_summary.update(fg)

        row = {
            "chain": chain,
            "resseq": resseq,
            "resname": resname,
            "fg_list": ";".join([k for k, v in fg.items() if v > 0]),
        }
        for k, v in fg.items():
            row[k] = int(v)
        annotated_rows.append(row)

    # Write annotated.csv
    extra_keys = sorted({k for r in annotated_rows for k in r.keys()
                         if k not in ("chain", "resseq", "resname", "fg_list") and k not in BASE_FG_KEYS})
    fg_keys = list(BASE_FG_KEYS) + extra_keys
    header = ["chain", "resseq", "resname", "fg_list"] + fg_keys

    with open(args.out_annotated, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in annotated_rows:
            full = {"chain": r.get("chain"), "resseq": r.get("resseq"), "resname": r.get("resname"), "fg_list": r.get("fg_list", "")}
            for k in fg_keys:
                full[k] = int(r.get(k, 0))
            w.writerow(full)

    # Write summary.csv
    with open(args.out_summary, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["functional_group", "count"])
        w.writeheader()
        for k in fg_keys:
            w.writerow({"functional_group": k, "count": int(all_summary.get(k, 0))})

    print(f"[OK] atoms (with H): {len(atoms)}")
    print(f"[OK] openbabel: {'ON' if args.use_openbabel and bond_order_map else 'OFF'}")
    print(f"[OK] residues annotated: {len(annotated_rows)}")
    print(f"[OK] wrote: {args.out_annotated}")
    print(f"[OK] wrote: {args.out_summary}")


if __name__ == "__main__":
    main()