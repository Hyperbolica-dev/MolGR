from __future__ import annotations

from pathlib import Path


PIPELINE_INIT = Path("src/molgr/_core/pipeline/__init__.pyi")


def _extract_class_block(lines: list[str], class_name: str) -> list[str]:
    start = None
    for i, line in enumerate(lines):
        if line.startswith(f"class {class_name}:"):
            start = i
            break
    if start is None:
        raise RuntimeError(f"missing class block: {class_name}")

    end = len(lines)
    for j in range(start + 1, len(lines)):
        if lines[j].startswith("class ") or lines[j].startswith("def "):
            end = j
            break
    return lines[start:end]


def _moduleize_class_block(block: list[str], class_name: str) -> str:
    body = block[1:]
    if body and body[0].strip() == '"""':
        doc_end = None
        for i in range(1, len(body)):
            if body[i].strip() == '"""':
                doc_end = i
                break
        if doc_end is not None:
            body = body[doc_end + 1 :]
    out: list[str] = []
    for raw in body:
        line = raw[4:] if raw.startswith("    ") else raw
        if line == "@staticmethod":
            continue
        line = line.replace(f"{class_name}.", "")
        out.append(line)
    return "\n".join(out).strip() + "\n"


def _extract_top_level_defs(lines: list[str]) -> str:
    start = None
    for i, line in enumerate(lines):
        if line.startswith("def "):
            start = i
            break
    if start is None:
        raise RuntimeError("missing top-level defs in pipeline stub")
    return "\n".join(lines[start:]).strip() + "\n"


def _write_module(path: Path, header: str, body: str) -> None:
    path.write_text(header + "\n" + body, encoding="utf-8")


def main() -> None:
    if not PIPELINE_INIT.exists():
        raise RuntimeError(f"missing pipeline init stub: {PIPELINE_INIT}")

    lines = PIPELINE_INIT.read_text(encoding="utf-8").splitlines()

    block_with_metals = _extract_class_block(lines, "reconstruct_with_metals")
    block_without_metals = _extract_class_block(lines, "reconstruct_without_metals")
    block_resonance = _extract_class_block(lines, "resonance")
    top_defs = _extract_top_level_defs(lines)

    with_metals_body = _moduleize_class_block(block_with_metals, "reconstruct_with_metals")
    without_metals_body = _moduleize_class_block(block_without_metals, "reconstruct_without_metals")
    resonance_body = _moduleize_class_block(block_resonance, "resonance")

    _write_module(
        Path("src/molgr/_core/pipeline/reconstruct_with_metals.pyi"),
        "\n".join(
            [
                '"""',
                "Fallback-aligned reconstruction helpers with metals",
                '"""',
                "",
                "from __future__ import annotations",
                "",
                "import collections.abc",
                "import typing",
                "",
                "import molgr._core.utils",
                "",
            ]
        ),
        with_metals_body,
    )

    _write_module(
        Path("src/molgr/_core/pipeline/reconstruct_without_metals.pyi"),
        "\n".join(
            [
                '"""',
                "Fallback-aligned no-metal reconstruction helpers",
                '"""',
                "",
                "from __future__ import annotations",
                "",
                "import typing",
                "",
                "import molgr._core.utils",
                "",
            ]
        ),
        without_metals_body,
    )

    _write_module(
        Path("src/molgr/_core/pipeline/resonance.pyi"),
        "\n".join(
            [
                '"""',
                "Fallback-aligned resonance helpers",
                '"""',
                "",
                "from __future__ import annotations",
                "",
                "import typing",
                "",
            ]
        ),
        resonance_body,
    )

    new_init = "\n".join(
        [
            '"""',
            "Pipeline-level helpers",
            '"""',
            "",
            "from __future__ import annotations",
            "",
            "import typing",
            "",
            "from . import reconstruct_with_metals, reconstruct_without_metals, resonance",
            "",
            "__all__: list[str] = [",
            '    "get_last_run_timing_breakdown_ms",',
            '    "get_radical_resonances_ptr",',
            '    "process_resonance_ptr",',
            '    "reconstruct_with_metals",',
            '    "reconstruct_without_metals",',
            '    "resonance",',
            '    "smiles_token_ptr",',
            "]",
            "",
            top_defs.strip(),
            "",
        ]
    )
    PIPELINE_INIT.write_text(new_init, encoding="utf-8")


if __name__ == "__main__":
    main()
