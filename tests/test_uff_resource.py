import os
from pathlib import Path

import molgr


def test_uff_parameter_file_is_packaged() -> None:
    resource = Path(molgr.__file__).resolve().parent / "data" / "UFF.prm"

    assert resource.is_file()
    assert resource.stat().st_size > 0


def test_pep561_marker_is_packaged() -> None:
    assert Path(molgr.__file__).resolve().with_name("py.typed").is_file()


def test_delvewheel_entrypoint_is_ascii() -> None:
    Path(molgr.__file__).read_bytes().decode("ascii")


def test_openbabel_data_dir_contains_uff_parameters() -> None:
    configured_dir = os.environ.get("BABEL_DATADIR")

    assert configured_dir is not None
    assert (Path(configured_dir) / "UFF.prm").is_file()
