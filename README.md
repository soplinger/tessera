# Tessera

**Offline wallpaper-group analyser: give it one small tile, get its symmetry group, its motif, and a picture of the repeat.**

Tessera reads a plain-text tile (a grid of characters that repeats forever in the plane), finds every rotation, mirror and glide that maps the pattern onto itself, and names its wallpaper group. It is useful for checking tile designs, quilt and mosaic layouts, textile repeats, and for teaching the 17 wallpaper groups. It runs on the Python standard library alone, with no servers, accounts, network calls or API keys.

## Features

- **Plain-text tile files** with `#` comments, `@name`, `@expect` and `@note` directives, and one character per cell.
- **Line-numbered lint**: ragged rows, stray whitespace, unknown directives, unknown group names and oversized tiles are all reported as `file:line: error: ...`.
- **Lattice analysis**: the translation lattice in reduced form, with its Bravais type (oblique, rectangular, centred or square) and the cell area.
- **Full symmetry enumeration** of the tile: every rotation, mirror and glide, with rotation centres, modulo translations.
- **Mirror versus glide** detection for each reflection direction, taking the tile's real translations into account.
- **Wallpaper group naming** for the 12 groups that square-pixel tiles can carry: p1, p2, pm, pg, cm, pmm, pmg, pgg, cmm, p4, p4m and p4g.
- **Motif finder**: the orbits of cells under the tile's symmetries, drawn as labels, so you can see the minimal block that makes the whole pattern.
- **Repeat warnings**: lint tells you when the pattern repeats in a smaller cell than the tile you wrote.
- **ASCII preview** of the tile repeated across and down.
- **SVG export** of the tiled pattern, with one colour per symbol.
- **Seeded generator**: `tessera generate pmg --seed 7` makes a random tile that classifies as exactly that group, retrying seeds until it does.
- **CI gate**: `tessera check` exits 0 on a match, 1 on a mismatch and 2 on a broken file, and reads the expectation from `@expect` or `--expect`.
- **JSON output** from `inspect --json` for scripts and dashboards.
- **Group reference**: `tessera groups` lists every group Tessera can name, with a plain description.

## Quick start

The repository ships a checkerboard example:

```text
$ python3 tessera.py inspect examples/checker.tile
examples/checker.tile: 2x2 tile 'Checkerboard'
group:    p4m (quarter-turns with mirrors through the rotation centres)
lattice:  square, basis (1,1) (-1,1)
period:   2 cells per repeat of 4
motif:    2 distinct cell orbits
symmetry: 8 point operations, 2 translations per tile
  ...one line per symmetry follows...
note:     Rotations about a cell centre and mirrors through cell centres keep the colours.
```

Then check it, preview it, and draw it:

```bash
python3 tessera.py check examples/checker.tile          # PASS ... p4m
python3 tessera.py motif examples/checker.tile           # two orbit labels: ab / ba
python3 tessera.py preview examples/checker.tile --across 4 --down 2
python3 tessera.py svg examples/checker.tile -o checker.svg
```

Make a tile from scratch and save it:

```bash
python3 tessera.py generate p4g --size 8 --seed 1 -o my-tile.tile
python3 tessera.py inspect my-tile.tile
```

Lint a broken file to see line-numbered errors:

```bash
python3 tessera.py lint examples/broken.tile
# examples/broken.tile:3: error: unknown directive '@colour'
# examples/broken.tile:6: error: row has 3 cells, expected 2 from line 4
```

## Tile file format

```text
# comments start with '#'
@name Optional display name
@expect p4m            # optional: the group a CI gate should require
@note free text, repeatable
ox                     # every other non-blank line is one row of the tile
xo
```

Rows must all have the same length, contain no spaces or tabs, and be at most 24 cells in each direction. Any printable, non-space character is a valid symbol.

## Commands

| Command | What it does |
| --- | --- |
| `groups` | List the wallpaper groups Tessera can name |
| `lint FILE` | Check a tile file; exit 2 on errors, print warnings otherwise |
| `inspect FILE [--json]` | Name the group and list its symmetries |
| `motif FILE` | Show the orbit label of every cell |
| `preview FILE [--across N] [--down N]` | Print the tile repeated as ASCII |
| `svg FILE -o OUT [--repeat N]` | Write an SVG picture of the tiled pattern |
| `generate GROUP [--size N] [--seed S] [--alphabet CHARS] [-o OUT]` | Make a random tile with exactly that group |
| `check FILE [--expect GROUP]` | CI gate: 0 on match, 1 on mismatch, 2 on errors |

## Scope and limits

- Tiles live on a square pixel grid, so only the 12 groups with 2-fold, 4-fold or mirror-only symmetry can appear. The hexagonal and triangular groups (p3, p3m1, p31m, p6, p6m) need 3-fold or 6-fold rotations, which cannot map square cells onto square cells. Tessera does not name them.
- Tiles are limited to 24 by 24 cells. Symmetry search cost grows with the square of the cell count, so keep tiles small for quick answers.
- `generate` builds the pattern by joining cells that the group's generators map onto each other, then picks random colours per joined class. It retries seeds until the output classifies as exactly the requested group.

## How to run

Requires Python 3.9 or newer. There is nothing to install.

```bash
python3 tessera.py --help
python3 tessera.py <command> ...
```

## How to test

Tests use pytest, which is the only third-party dependency.

```bash
python3 -m pip install pytest
python3 -m pytest -q
```

The tests cover lattice arithmetic, the tile parser's line-numbered errors, the symmetry and group classification of known patterns, a round trip for every generated group, rendering, and the CLI exit codes.

## Status

The code and tests were written in an authoring session that could not execute Python, so the test suite has not yet been run. Running `python3 -m pytest -q` is the first check to do.

## Licence

No licence file is included; add one before publishing.
