"""
Author: TMJ
Date: 2026-02-25 00:49:31
LastEditors: TMJ
LastEditTime: 2026-02-25 13:45:37
Description: 请填写简介
"""
# pyright: reportMissingImports=false

import pytest


pytest.importorskip("openbabel")

from molgr.fallback.utils.consts import get_possible_metal_radicals


def test_get_possible_metal_radicals_unknown_metal_returns_empty_set() -> None:
    assert get_possible_metal_radicals("UUnobtainium", 1) == set()


def test_get_possible_metal_radicals_high_valence_returns_empty_set() -> None:
    assert get_possible_metal_radicals("Li", 10) == set()


def test_get_possible_metal_radicals_known_cpp_parity_cases() -> None:
    assert get_possible_metal_radicals("Fe", 2) == {4, 2, 0}
    assert get_possible_metal_radicals("Cu", 2) == {1}


def test_get_possible_metal_radicals_allows_d_shell_rearrangement_for_co_i() -> None:
    assert get_possible_metal_radicals("Co", 1) == {2, 0}


def test_get_possible_metal_radicals_f_shell_branch_tracks_remaining_parity() -> None:
    assert get_possible_metal_radicals("Ce", 4) == {0}
