import json

import pytest

import tessera
from tessera import Tile, analyse, find_generated, parse_tile


def _det(basis):
    (a, b), (c, d) = basis
    return abs(a * d - b * c)


def test_d4_has_eight_isometries():
    assert len(tessera.D4) == 8


def test_lattice_basis_of_full_period():
    basis = tessera.lattice_basis([(0, 0), (8, 0), (0, 8)])
    assert basis == ((8, 0), (0, 8))


def test_lattice_basis_picks_up_half_shift():
    # (4, 4) plus the 8x8 period spans a lattice of determinant 32.
    basis = tessera.lattice_basis([(0, 0), (4, 4), (8, 0), (0, 8)])
    assert _det(basis) == 32


def test_reduced_centred_pair_is_square():
    u, v = tessera.reduce_basis((4, 4), (0, 8))
    assert {u, v} == {(4, 4), (-4, 4)}
    assert tessera.lattice_kind(u, v)[0] == "square"


def test_rectangular_lattice_is_rect():
    assert tessera.lattice_kind((8, 0), (0, 4))[0] == "rect"


def test_skewed_lattice_is_oblique():
    assert tessera.lattice_kind((1, 2), (2, -2))[0] == "oblique"


def test_parse_skips_comments_and_blank_lines():
    tile = parse_tile("# hello\n\n@name Demo\nab\nba\n")
    assert tile.errors == []
    assert tile.rows == ["ab", "ba"]
    assert tile.name == "Demo"


def test_parse_reports_ragged_row_with_line_number():
    tile = parse_tile("ab\nabc\n")
    assert len(tile.errors) == 1
    assert tile.errors[0].line == 2
    assert "row has 3 cells" in tile.errors[0].message


def test_parse_reports_unknown_directive():
    tile = parse_tile("@colour red\nab\n")
    assert [(i.line, "unknown directive" in i.message) for i in tile.errors] == [(1, True)]


def test_parse_rejects_unknown_expected_group():
    tile = parse_tile("@expect p3\nab\n")
    assert tile.errors[0].line == 1
    assert "unknown wallpaper group 'p3'" in tile.errors[0].message


def test_parse_rejects_whitespace_inside_row():
    tile = parse_tile("a b\n")
    assert tile.errors[0].line == 1
    assert "column 2" in tile.errors[0].message


def test_parse_rejects_file_without_rows():
    tile = parse_tile("# only a comment\n")
    assert tile.errors[0].line == 0
    assert "no tile rows" in tile.errors[0].message


def test_uniform_tile_has_full_symmetry():
    assert analyse(Tile(rows=["oo", "oo"])).group == "p4m"


def test_uniform_tile_is_one_motif_orbit():
    assert analyse(Tile(rows=["oo", "oo"])).motif_count == 1


def test_checkerboard_is_p4m_with_two_orbits():
    an = analyse(Tile(rows=["ox", "xo"]))
    assert an.group == "p4m"
    assert an.period_cells == 2
    assert an.motif_count == 2


@pytest.mark.parametrize("group", list(tessera.GROUPS))
def test_generated_tile_classifies_as_requested_group(group):
    rows = find_generated(group, 8, 1, "ox+*")
    assert analyse(Tile(rows=rows)).group == group


def test_generated_pg_has_a_glide_line():
    rows = find_generated("pg", 8, 1, "ox+*")
    lines = tessera.describe_operations(analyse(Tile(rows=rows)))
    assert any(line.startswith("glide along (1, 0)") for line in lines)


def test_render_ascii_repeats_rows_and_columns():
    out = tessera.render_ascii(["ab", "cd"], 2, 2)
    assert out == ["abab", "cdcd", "abab", "cdcd"]


def test_render_svg_draws_every_repeated_cell():
    svg = tessera.render_svg(["ab", "ba"], repeat=2)
    assert svg.count("<rect") == 2 * 2 * 4
    assert svg.startswith("<svg")


def test_generate_command_rejects_odd_size(capsys):
    assert tessera.main(["generate", "p1", "--size", "7"]) == 2
    assert "--size must be an even number" in capsys.readouterr().err


def test_lint_warns_when_tile_repeats_early(tmp_path, capsys):
    path = tmp_path / "rep.tile"
    path.write_text("abab\nabab\n", encoding="utf-8")
    assert tessera.main(["lint", str(path)]) == 0
    assert "smaller tile" in capsys.readouterr().out


def test_check_passes_and_fails_with_exit_codes(tmp_path, capsys):
    path = tmp_path / "ok.tile"
    path.write_text("@expect p4m\noo\noo\n", encoding="utf-8")
    assert tessera.main(["check", str(path)]) == 0
    assert tessera.main(["check", str(path), "--expect", "pm"]) == 1
    assert "FAIL" in capsys.readouterr().out


def test_check_needs_an_expectation(tmp_path):
    path = tmp_path / "none.tile"
    path.write_text("oo\noo\n", encoding="utf-8")
    assert tessera.main(["check", str(path)]) == 2


def test_check_reports_parse_errors_with_exit_two(tmp_path, capsys):
    path = tmp_path / "bad.tile"
    path.write_text("ab\nabc\n", encoding="utf-8")
    assert tessera.main(["check", str(path), "--expect", "p1"]) == 2
    assert f"{path}:2: error:" in capsys.readouterr().out


def test_generate_command_emits_matching_expect_line(capsys):
    assert tessera.main(["generate", "pmm", "--size", "8", "--seed", "3"]) == 0
    text = capsys.readouterr().out
    assert text.startswith("@name generated pmm\n@expect pmm\n")
    assert analyse(parse_tile(text)).group == "pmm"


def test_inspect_json_names_the_group(tmp_path, capsys):
    path = tmp_path / "p4m.tile"
    path.write_text("oo\noo\n", encoding="utf-8")
    assert tessera.main(["inspect", "--json", str(path)]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["group"] == "p4m"
    assert payload["tile_cells"] == 4
