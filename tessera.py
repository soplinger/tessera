#!/usr/bin/env python3
"""Tessera: an offline wallpaper-group analyser for square-pixel tile files.

A tile file is a small grid of single characters that repeats forever in the
plane. Tessera finds every rotation, mirror and glide that maps the pattern
onto itself, names its wallpaper group, and can draw, generate and gate it.
Standard library only; nothing touches the network.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product

MAX_SIZE = 24

Vec = tuple[int, int]
Mat = tuple[Vec, Vec]

IDENTITY: Mat = ((1, 0), (0, 1))
NEG_I: Mat = ((-1, 0), (0, -1))
ROT90: Mat = ((0, -1), (1, 0))
ROT270: Mat = ((0, 1), (-1, 0))
MIRROR_X: Mat = ((1, 0), (0, -1))  # y -> -y, mirror axis along (1, 0)
MIRROR_Y: Mat = ((-1, 0), (0, 1))  # x -> -x, mirror axis along (0, 1)
MIRROR_D: Mat = ((0, 1), (1, 0))  # swap x and y, axis along (1, 1)
MIRROR_AD: Mat = ((0, -1), (-1, 0))  # axis along (1, -1)

# Mirror axis direction of each reflection. The same vector is the functional
# whose kernel is the mirror's normal direction, which is what the mirror test uses.
AXIS: dict[Mat, Vec] = {
    MIRROR_X: (1, 0),
    MIRROR_Y: (0, 1),
    MIRROR_D: (1, 1),
    MIRROR_AD: (1, -1),
}

# The eight integer isometries of the square grid (the point group D4).
D4: list[Mat] = sorted(
    m
    for m in (
        ((a, b), (c, d)) for a, b, c, d in product((-1, 0, 1), repeat=4)
    )
    if m[0][0] ** 2 + m[0][1] ** 2 == 1
    and m[1][0] ** 2 + m[1][1] ** 2 == 1
    and m[0][0] * m[1][0] + m[0][1] * m[1][1] == 0
)

GROUPS: dict[str, str] = {
    "p1": "translations only; no rotations or reflections",
    "p2": "half-turn rotations only, no mirrors or glides",
    "pm": "parallel mirror lines on a primitive rectangular cell",
    "pg": "parallel glide lines only, on a primitive rectangular cell",
    "cm": "mirror lines on a centred rectangular cell",
    "pmm": "mirrors in two perpendicular directions on a rectangular cell",
    "pmg": "mirrors in one direction, perpendicular glides, rectangular cell",
    "pgg": "perpendicular glide lines only, rectangular cell",
    "cmm": "mirrors in two perpendicular directions on a centred cell",
    "p4": "quarter-turn rotations on a square cell, no reflections",
    "p4m": "quarter-turns with mirrors through the rotation centres",
    "p4g": "quarter-turns, diagonal mirrors and side glides on a square cell",
}

# Canonical generators on an N x N cell (N even, h = N // 2). Each entry is
# (linear part, translation). Used by `generate` and by the round-trip tests.
def _generators(h: int) -> dict[str, list[tuple[Mat, Vec]]]:
    return {
        "p1": [],
        "p2": [(NEG_I, (0, 0))],
        "pm": [(MIRROR_X, (0, 0))],
        "pg": [(MIRROR_X, (h, 0))],
        "cm": [(MIRROR_X, (0, 0)), (IDENTITY, (h, h))],
        "pmm": [(MIRROR_X, (0, 0)), (MIRROR_Y, (0, 0))],
        "pmg": [(MIRROR_X, (0, 0)), (MIRROR_Y, (0, h))],
        "pgg": [(MIRROR_X, (h, 0)), (MIRROR_Y, (0, h))],
        "cmm": [(MIRROR_X, (0, 0)), (MIRROR_Y, (0, 0)), (IDENTITY, (h, h))],
        "p4": [(ROT90, (0, 0))],
        "p4m": [(ROT90, (0, 0)), (MIRROR_X, (0, 0))],
        "p4g": [(ROT90, (0, 0)), (MIRROR_D, (h, -h))],
    }


# ---------------------------------------------------------------- vectors

def _apply(a: Mat, p: Vec) -> Vec:
    return (a[0][0] * p[0] + a[0][1] * p[1], a[1][0] * p[0] + a[1][1] * p[1])


def _dot(a: Vec, b: Vec) -> int:
    return a[0] * b[0] + a[1] * b[1]


def _cross(a: Vec, b: Vec) -> int:
    return a[0] * b[1] - a[1] * b[0]


def _egcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, s, t) with s*a + t*b == g == gcd(a, b), for a, b >= 0."""
    if b == 0:
        return a, 1, 0
    g, s, t = _egcd(b, a % b)
    return g, t, s - (a // b) * t


def lattice_basis(vectors: list[Vec]) -> tuple[Vec, Vec]:
    """Hermite-style basis of the integer lattice spanned by `vectors`.

    The vectors must span a full-rank lattice (the tile always adds (W, 0)
    and (0, H), so they do). Returns ((a, b), (0, d)) with a > 0, d > 0.
    """
    x_vec: Vec | None = None
    zeros: list[Vec] = []
    for v in vectors:
        if v[0] == 0:
            zeros.append(v)
            continue
        if x_vec is None:
            x_vec = v
            continue
        a, b = x_vec[0], v[0]
        g, s, t = _egcd(abs(a), abs(b))
        s = s if a > 0 else -s
        t = t if b > 0 else -t
        new = (s * x_vec[0] + t * v[0], s * x_vec[1] + t * v[1])
        other = (
            (b // g) * x_vec[0] - (a // g) * v[0],
            (b // g) * x_vec[1] - (a // g) * v[1],
        )
        zeros.append(other)
        x_vec = new
    assert x_vec is not None
    if x_vec[0] < 0:
        x_vec = (-x_vec[0], -x_vec[1])
    d = 0
    for z in zeros:
        d = math.gcd(d, abs(z[1]))
    return (x_vec[0], x_vec[1] % d), (0, d)


def reduce_basis(u: Vec, v: Vec) -> tuple[Vec, Vec]:
    """Gauss-Lagrange reduction: shortest-first basis with |u·v| <= |u|²/2."""
    while True:
        if _dot(u, u) > _dot(v, v):
            u, v = v, u
        nu = _dot(u, u)
        m = (2 * _dot(u, v) + nu) // (2 * nu)  # round(u·v / |u|²)
        if m == 0:
            return u, v
        v = (v[0] - m * u[0], v[1] - m * u[1])


def lattice_kind(u: Vec, v: Vec) -> tuple[str, tuple[Vec, Vec] | None]:
    """Bravais type of a reduced 2-D lattice, plus its rectangular axes.

    'square' and 'rect' have an orthogonal primitive basis. 'centred' has an
    orthogonal pair of lattice vectors spanning half the cell. Anything else is
    'oblique'.
    """
    det = abs(_cross(u, v))
    if _dot(u, v) == 0:
        kind = "square" if _dot(u, u) == _dot(v, v) else "rect"
        return kind, (u, v)
    bound = 2 * max(abs(c) for c in (*u, *v)) + 2
    found = [
        (x, y)
        for x in range(-bound, bound + 1)
        for y in range(-bound, bound + 1)
        if (x, y) != (0, 0)
        and _cross((x, y), v) % det == 0
        and _cross(u, (x, y)) % det == 0
    ]
    for i, a in enumerate(found):
        for b in found[i + 1:]:
            if _dot(a, b) == 0 and abs(_cross(a, b)) == 2 * det:
                return "centred", (a, b)
    return "oblique", None


def _aligned(a: Mat, axes: tuple[Vec, Vec] | None) -> bool:
    """True when the reflection's mirror axis is parallel to a lattice axis."""
    if axes is None:
        return False
    axis = AXIS[a]
    return any(_cross(axis, e) == 0 for e in axes)


# ---------------------------------------------------------------- tile files

@dataclass
class Issue:
    line: int
    message: str


@dataclass
class Tile:
    rows: list[str]
    name: str = ""
    expect: str = ""
    notes: list[str] = field(default_factory=list)
    errors: list[Issue] = field(default_factory=list)

    @property
    def width(self) -> int:
        return len(self.rows[0])

    @property
    def height(self) -> int:
        return len(self.rows)

    def cell(self, x: int, y: int) -> str:
        return self.rows[y % self.height][x % self.width]


def parse_tile(text: str) -> Tile:
    """Parse a tile file. Errors carry line numbers and never raise."""
    rows: list[str] = []
    row_lines: list[int] = []
    errors: list[Issue] = []
    name, expect, expect_line = "", "", 0
    notes: list[str] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("@"):
            key, _, value = stripped[1:].partition(" ")
            value = value.strip()
            if key == "name":
                name = value
            elif key == "expect":
                expect, expect_line = value, lineno
            elif key == "note":
                notes.append(value)
            else:
                errors.append(Issue(lineno, f"unknown directive '@{key}'"))
            continue
        blank = next((i for i, ch in enumerate(raw) if ch.isspace()), None)
        if blank is not None:
            errors.append(
                Issue(lineno, f"whitespace in tile row at column {blank + 1}")
            )
            continue
        rows.append(stripped)
        row_lines.append(lineno)

    if not rows:
        errors.append(Issue(0, "no tile rows found"))
        return Tile(rows=[" "], name=name, expect=expect, notes=notes, errors=errors)

    width = len(rows[0])
    for line, row in zip(row_lines, rows):
        if len(row) != width:
            errors.append(
                Issue(line, f"row has {len(row)} cells, expected {width} from line {row_lines[0]}")
            )
    if width > MAX_SIZE or len(rows) > MAX_SIZE:
        errors.append(
            Issue(row_lines[0], f"tile is {width}x{len(rows)}; the largest supported is {MAX_SIZE}x{MAX_SIZE}")
        )
    if expect and expect not in GROUPS:
        errors.append(
            Issue(expect_line, f"unknown wallpaper group '{expect}' (see `tessera groups`)")
        )
    if errors:
        return Tile(rows=rows, name=name, expect=expect, notes=notes, errors=errors)
    return Tile(rows=rows, name=name, expect=expect, notes=notes)


# ---------------------------------------------------------------- symmetry

def _in_period(p: Vec, w: int, h: int) -> bool:
    return p[0] % w == 0 and p[1] % h == 0


def _symmetric(tile: Tile, a: Mat, t: Vec) -> bool:
    """Does the map i -> a·i + t leave the periodic pattern unchanged?"""
    w, h = tile.width, tile.height
    if _in_period(_apply(a, (w, 0)), w, h) and _in_period(_apply(a, (0, h)), w, h):
        kx, ky = w, h
    else:
        # The map does not preserve the tile's period lattice, so the check
        # must cover lcm(w, h) cells to see every residue.
        kx = ky = w * h // math.gcd(w, h)
    for x in range(kx):
        for y in range(ky):
            px, py = _apply(a, (x, y))
            if tile.cell(px + t[0], py + t[1]) != tile.cell(x, y):
                return False
    return True


def _find(parent: list[int], i: int) -> int:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def _union(parent: list[int], i: int, j: int) -> None:
    ri, rj = _find(parent, i), _find(parent, j)
    if ri != rj:
        parent[max(ri, rj)] = min(ri, rj)


_ORBIT_LABELS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"


@dataclass
class Analysis:
    group: str | None
    lattice: str
    basis: tuple[Vec, Vec]
    period_cells: int
    tile_cells: int
    elements: list[tuple[Mat, Vec]]
    mirrors: dict[Mat, bool]
    motif: list[list[int]]
    motif_count: int

    @property
    def point_group(self) -> list[Mat]:
        return [a for a in D4 if any(b == a for b, _ in self.elements)]

    @property
    def translations(self) -> int:
        return sum(1 for a, _ in self.elements if a == IDENTITY)


def _classify(
    point_group: list[Mat],
    mirrors: dict[Mat, bool],
    lattice: str,
    axes: tuple[Vec, Vec] | None,
) -> str | None:
    if ROT90 in point_group or ROT270 in point_group:
        reflections = [a for a in point_group if a in mirrors]
        if not reflections:
            return "p4"
        sides = [a for a in reflections if _aligned(a, axes)]
        return "p4m" if any(mirrors[a] for a in sides) else "p4g"
    reflections = [a for a in point_group if a in mirrors]
    if not reflections:
        return "p2" if NEG_I in point_group else "p1"
    # A square lattice whose mirrors all sit on the diagonals is a centred
    # rectangle once the mirrors are taken as the axes.
    centred = lattice == "centred" or (
        lattice == "square" and not any(_aligned(a, axes) for a in reflections)
    )
    if len(reflections) == 1:
        a = reflections[0]
        if centred:
            return "cm" if mirrors[a] else None
        return "pm" if mirrors[a] else "pg"
    if len(reflections) == 2:
        if centred:
            return "cmm" if any(mirrors.values()) else None
        return {2: "pmm", 1: "pmg", 0: "pgg"}[sum(mirrors[a] for a in reflections)]
    return None


def analyse(tile: Tile) -> Analysis:
    """Find every symmetry of the tile's pattern and name its wallpaper group."""
    w, h = tile.width, tile.height
    translations = [
        (x, y) for x in range(w) for y in range(h) if _symmetric(tile, IDENTITY, (x, y))
    ]
    basis = lattice_basis(translations + [(w, 0), (0, h)])
    u, v = reduce_basis(*basis)
    lattice, axes = lattice_kind(u, v)
    elements = [
        (a, (x, y))
        for a in D4
        for x in range(w)
        for y in range(h)
        if _symmetric(tile, a, (x, y))
    ]
    point = [a for a in D4 if any(b == a for b, _ in elements)]
    mirrors: dict[Mat, bool] = {}
    for a in point:
        if a in AXIS:
            t0 = next(t for b, t in elements if b == a)
            g = math.gcd(abs(_dot(AXIS[a], u)), abs(_dot(AXIS[a], v)))
            # A mirror exists when some lattice shift cancels the glide part.
            mirrors[a] = _dot(AXIS[a], t0) % g == 0
    motif, count = _motif(tile, elements)
    period = abs(_cross(u, v))
    return Analysis(
        group=_classify(point, mirrors, lattice, axes),
        lattice=lattice,
        basis=(u, v),
        period_cells=period,
        tile_cells=w * h,
        elements=elements,
        mirrors=mirrors,
        motif=motif,
        motif_count=count,
    )


def _motif(tile: Tile, elements: list[tuple[Mat, Vec]]) -> tuple[list[list[int]], int]:
    """Label each cell by its orbit under the tile-preserving symmetries."""
    w, h = tile.width, tile.height
    parent = list(range(w * h))
    for a, t in elements:
        if not (_in_period(_apply(a, (w, 0)), w, h) and _in_period(_apply(a, (0, h)), w, h)):
            continue
        for x in range(w):
            for y in range(h):
                px, py = _apply(a, (x, y))
                px, py = (px + t[0]) % w, (py + t[1]) % h
                _union(parent, y * w + x, py * w + px)
    labels: dict[int, int] = {}
    grid: list[list[int]] = []
    for y in range(h):
        row = []
        for x in range(w):
            root = _find(parent, y * w + x)
            if root not in labels:
                labels[root] = len(labels)
            row.append(labels[root])
        grid.append(row)
    return grid, len(labels)


def _centre(a: Mat, t: Vec) -> tuple[Fraction, Fraction] | None:
    """Fixed point of x -> a·x + t, if one exists."""
    m00, m01 = 1 - a[0][0], -a[0][1]
    m10, m11 = -a[1][0], 1 - a[1][1]
    det = m00 * m11 - m01 * m10
    if det == 0:
        return None
    return (
        Fraction(m11 * t[0] - m01 * t[1], det),
        Fraction(-m10 * t[0] + m00 * t[1], det),
    )


def _op_name(a: Mat) -> str:
    if a == IDENTITY:
        return "translation"
    if a == NEG_I:
        return "2-fold rotation"
    if a in (ROT90, ROT270):
        return "4-fold rotation"
    return "reflection"


def describe_operations(an: Analysis) -> list[str]:
    """One line per non-translation symmetry, positions taken mod the tile."""
    lines = []
    for a, t in an.elements:
        if a == IDENTITY:
            continue
        if a in AXIS:
            kind = "mirror" if an.mirrors[a] else "glide"
            axis = AXIS[a]
            lines.append(f"{kind} along ({axis[0]}, {axis[1]}), shift {t[0]},{t[1]}")
        else:
            c = _centre(a, t)
            where = f"({c[0]}, {c[1]})" if c else "?"
            lines.append(f"{_op_name(a)} about {where}")
    return lines


# ---------------------------------------------------------------- lint

def lint(tile: Tile) -> list[str]:
    """Warnings about a valid tile. Errors are already on tile.errors."""
    warnings = []
    if len(set("".join(tile.rows))) == 1:
        warnings.append("every cell is the same symbol; the pattern has maximal symmetry")
    an = analyse(tile)
    if an.period_cells < an.tile_cells:
        warnings.append(
            f"the pattern repeats every {an.period_cells} cells but the tile holds "
            f"{an.tile_cells}; a smaller tile would describe it"
        )
    return warnings


# ---------------------------------------------------------------- generation

def generate(group: str, size: int, seed: int, alphabet: str) -> list[str]:
    """Random tile with exactly the requested symmetry, from one seed."""
    n = size
    parent = list(range(n * n))
    for a, t in _generators(n // 2)[group]:
        for x in range(n):
            for y in range(n):
                px, py = _apply(a, (x, y))
                px, py = (px + t[0]) % n, (py + t[1]) % n
                _union(parent, y * n + x, py * n + px)
    rng = random.Random(seed)
    roots = sorted({_find(parent, i) for i in range(n * n)})
    colour = {r: rng.choice(alphabet) for r in roots}
    return [
        "".join(colour[_find(parent, y * n + x)] for x in range(n))
        for y in range(n)
    ]


def find_generated(group: str, size: int, seed: int, alphabet: str) -> list[str]:
    """Try seeds until the generated tile classifies as exactly `group`."""
    for s in range(seed, seed + 200):
        rows = generate(group, size, s, alphabet)
        if analyse(Tile(rows=rows)).group == group:
            return rows
    raise RuntimeError(f"could not generate a tile for {group} from seed {seed}")


# ---------------------------------------------------------------- rendering

def render_ascii(rows: list[str], across: int, down: int) -> list[str]:
    out = []
    for _ in range(down):
        for row in rows:
            out.append(row * across)
    return out


_PALETTE = ["#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#8ab17d", "#b56576", "#6d6875"]


def render_svg(rows: list[str], repeat: int = 3, cell: int = 16) -> str:
    chars = sorted(set("".join(rows)))
    colour = {c: _PALETTE[i % len(_PALETTE)] for i, c in enumerate(chars)}
    w, h = len(rows[0]), len(rows)
    width, height = w * repeat * cell, h * repeat * cell
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    ]
    for ry in range(repeat):
        for rx in range(repeat):
            for y, row in enumerate(rows):
                for x, ch in enumerate(row):
                    px = (rx * w + x) * cell
                    py = (ry * h + y) * cell
                    parts.append(
                        f'<rect x="{px}" y="{py}" width="{cell}" height="{cell}" '
                        f'fill="{colour[ch]}"/>'
                    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


# ---------------------------------------------------------------- CLI

def _load(path: str) -> Tile | None:
    with open(path, encoding="utf-8") as fh:
        tile = parse_tile(fh.read())
    if tile.errors:
        for issue in tile.errors:
            if issue.line:
                print(f"{path}:{issue.line}: {issue.message}", file=sys.stderr)
            else:
                print(f"{path}: {issue.message}", file=sys.stderr)
        return None
    return tile


def _inspect_text(path: str, tile: Tile, an: Analysis) -> list[str]:
    u, v = an.basis
    lines = [
        f"{path}: {tile.width}x{tile.height} tile" + (f" '{tile.name}'" if tile.name else ""),
        f"group:    {an.group or 'unrecognised'}"
        + (f" ({GROUPS[an.group]})" if an.group in GROUPS else ""),
        f"lattice:  {an.lattice}, basis ({u[0]},{u[1]}) ({v[0]},{v[1]})",
        f"period:   {an.period_cells} cells per repeat of {an.tile_cells}",
        f"motif:    {an.motif_count} distinct cell orbits",
        f"symmetry: {len(an.point_group)} point operations, {an.translations} translations per tile",
    ]
    for line in describe_operations(an):
        lines.append(f"  {line}")
    for note in tile.notes:
        lines.append(f"note:     {note}")
    return lines


def cmd_groups(_: argparse.Namespace) -> int:
    for name, text in GROUPS.items():
        print(f"{name:4} {text}")
    print()
    print("Square-pixel tiles cannot carry 3-fold or 6-fold symmetry (p3, p3m1, p31m,")
    print("p6, p6m), so those groups are not reported.")
    return 0


def cmd_lint(args: argparse.Namespace) -> int:
    with open(args.file, encoding="utf-8") as fh:
        tile = parse_tile(fh.read())
    if tile.errors:
        for issue in tile.errors:
            loc = f"{args.file}:{issue.line}" if issue.line else args.file
            print(f"{loc}: error: {issue.message}")
        return 2
    for w in lint(tile):
        print(f"{args.file}: warning: {w}")
    print(f"{args.file}: ok ({tile.width}x{tile.height})")
    return 0


def cmd_inspect(args: argparse.Namespace) -> int:
    tile = _load(args.file)
    if tile is None:
        return 2
    an = analyse(tile)
    if args.json:
        payload = {
            "name": tile.name,
            "group": an.group,
            "lattice": an.lattice,
            "basis": [list(an.basis[0]), list(an.basis[1])],
            "period_cells": an.period_cells,
            "tile_cells": an.tile_cells,
            "motif_cells": an.motif_count,
            "translations": an.translations,
            "point_operations": len(an.point_group),
            "operations": describe_operations(an),
        }
        print(json.dumps(payload, indent=2))
    else:
        print("\n".join(_inspect_text(args.file, tile, an)))
    return 0


def cmd_motif(args: argparse.Namespace) -> int:
    tile = _load(args.file)
    if tile is None:
        return 2
    an = analyse(tile)
    print(f"{an.motif_count} orbit(s) in a {tile.width}x{tile.height} tile")
    for row in an.motif:
        print("".join(_ORBIT_LABELS[i] if i < len(_ORBIT_LABELS) else "*" for i in row))
    return 0


def cmd_preview(args: argparse.Namespace) -> int:
    tile = _load(args.file)
    if tile is None:
        return 2
    for line in render_ascii(tile.rows, args.across, args.down):
        print(line)
    return 0


def cmd_svg(args: argparse.Namespace) -> int:
    tile = _load(args.file)
    if tile is None:
        return 2
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(render_svg(tile.rows, repeat=args.repeat))
    print(f"wrote {args.output}")
    return 0


def cmd_generate(args: argparse.Namespace) -> int:
    if args.group not in GROUPS:
        print(f"unknown group '{args.group}' (see `tessera groups`)", file=sys.stderr)
        return 2
    if args.size % 2 or not 4 <= args.size <= MAX_SIZE:
        print(f"--size must be an even number from 4 to {MAX_SIZE}", file=sys.stderr)
        return 2
    if len(set(args.alphabet)) < 2 or any(c.isspace() for c in args.alphabet):
        print("--alphabet needs at least two distinct non-space characters", file=sys.stderr)
        return 2
    try:
        rows = find_generated(args.group, args.size, args.seed, args.alphabet)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    text = f"@name generated {args.group}\n@expect {args.group}\n" + "\n".join(rows) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {args.output}")
    else:
        sys.stdout.write(text)
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    with open(args.file, encoding="utf-8") as fh:
        tile = parse_tile(fh.read())
    if tile.errors:
        for issue in tile.errors:
            loc = f"{args.file}:{issue.line}" if issue.line else args.file
            print(f"{loc}: error: {issue.message}")
        return 2
    expected = args.expect or tile.expect
    if not expected:
        print(f"{args.file}: no expectation; pass --expect or add an @expect line", file=sys.stderr)
        return 2
    found = analyse(tile).group or "unrecognised"
    if found == expected:
        print(f"PASS {args.file}: {found}")
        return 0
    print(f"FAIL {args.file}: expected {expected}, found {found}")
    return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tessera",
        description="Offline wallpaper-group analyser for square-pixel tile files.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("groups", help="list the wallpaper groups Tessera can name")

    p = sub.add_parser("lint", help="check a tile file with line-numbered errors")
    p.add_argument("file")
    p.set_defaults(func=cmd_lint)

    p = sub.add_parser("inspect", help="name the group and list its operations")
    p.add_argument("file")
    p.add_argument("--json", action="store_true", help="print JSON instead of text")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("motif", help="show the minimal motif as orbit labels")
    p.add_argument("file")
    p.set_defaults(func=cmd_motif)

    p = sub.add_parser("preview", help="print the tile repeated as ASCII")
    p.add_argument("file")
    p.add_argument("--across", type=int, default=3)
    p.add_argument("--down", type=int, default=3)
    p.set_defaults(func=cmd_preview)

    p = sub.add_parser("svg", help="write an SVG picture of the tiled pattern")
    p.add_argument("file")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--repeat", type=int, default=3)
    p.set_defaults(func=cmd_svg)

    p = sub.add_parser("generate", help="make a random tile with a given group")
    p.add_argument("group")
    p.add_argument("--size", type=int, default=8)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--alphabet", default="ox+*")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_generate)

    p = sub.add_parser("check", help="CI gate: exit 0 only if the group matches")
    p.add_argument("file")
    p.add_argument("--expect", help="group to require (defaults to @expect)")
    p.set_defaults(func=cmd_check)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "groups":
        return cmd_groups(args)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
