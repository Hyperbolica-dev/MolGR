from pathlib import Path

import molgr


def test_uff_parameter_file_is_packaged() -> None:
    resource = Path(molgr.__file__).resolve().parent / "data" / "UFF.prm"

    assert resource.is_file()
    assert resource.stat().st_size > 0
