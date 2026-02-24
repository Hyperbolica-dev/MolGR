from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
from typing import Any, Optional
from urllib.parse import quote

from jinja2 import Environment, FileSystemLoader, select_autoescape


TRACE_FILENAME = "trace.jsonl"
REPORT_FILENAME = "report.html"
_TEMPLATE_NAME = "molgr_debug_trace_report.html.j2"
_TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"


_METRIC_NOOP_KEYS = ("charge", "radical", "atom_count", "bond_count")
_ALWAYS_VISIBLE_META_KEYS = {
    "reason",
    "cap",
    "cap_value",
    "truncated_by_cap",
    "candidate_count",
    "candidate_smiles",
    "traceback",
}


def _read_trace_events(
    trace_path: Path, *, max_events_rendered: int
) -> tuple[list[dict[str, Any]], int, int, bool]:
    events: list[dict[str, Any]] = []
    valid_events = 0
    malformed_lines = 0

    if not trace_path.exists():
        raise FileNotFoundError(f"trace file not found: {trace_path}")

    with trace_path.open("r", encoding="utf-8") as trace_fh:
        for line in trace_fh:
            raw = line.strip()
            if not raw:
                continue

            try:
                event = json.loads(raw)
            except json.JSONDecodeError:
                malformed_lines += 1
                continue

            if not isinstance(event, dict):
                malformed_lines += 1
                continue

            valid_events += 1
            if len(events) < max_events_rendered:
                events.append(event)

    truncated = valid_events > max_events_rendered
    return events, valid_events, malformed_lines, truncated


class _SmilesRenderer:
    def __init__(self, *, max_atoms_svg: int) -> None:
        self.max_atoms_svg = max_atoms_svg
        self._cache: dict[str, Optional[str]] = {}

    def render(self, smiles: Optional[str]) -> Optional[str]:
        if not isinstance(smiles, str) or smiles == "":
            return None
        if smiles in self._cache:
            return self._cache[smiles]

        svg: Optional[str] = None
        try:
            from openbabel import pybel  # type: ignore[import-not-found]

            pybel.ob.obErrorLog.StopLogging()

            m = pybel.readstring("smi", smiles)
            svg = None if m.OBMol.NumAtoms() > self.max_atoms_svg else m.write("svg")
        except Exception:
            svg = None

        self._cache[smiles] = svg
        return svg


def _meta_dict(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return dict(value)
    return {}


def _merge_meta(start_meta: Any, end_meta: Any) -> Optional[dict[str, Any]]:
    merged = _meta_dict(start_meta)
    merged.update(_meta_dict(end_meta))
    if not merged:
        return None
    return merged


def _same_metrics(metrics_before: Any, metrics_after: Any) -> bool:
    if not isinstance(metrics_before, dict) or not isinstance(metrics_after, dict):
        return True
    return all(metrics_before.get(key) == metrics_after.get(key) for key in _METRIC_NOOP_KEYS)


def _metric_delta_table(metrics_before: Any, metrics_after: Any) -> Optional[list[dict[str, Any]]]:
    if not isinstance(metrics_before, dict) or not isinstance(metrics_after, dict):
        return None

    keys: list[str] = []
    for key in _METRIC_NOOP_KEYS:
        if key in metrics_before or key in metrics_after:
            keys.append(key)
    for key in sorted(set(metrics_before.keys()) | set(metrics_after.keys())):
        if key in _METRIC_NOOP_KEYS:
            continue
        keys.append(str(key))

    out: list[dict[str, Any]] = []
    for key in keys:
        before = metrics_before.get(key)
        after = metrics_after.get(key)
        if before == after:
            continue
        delta: Optional[float] = None
        if isinstance(before, (int, float)) and isinstance(after, (int, float)):
            delta = float(after - before)
        out.append({"metric": key, "before": before, "after": after, "delta": delta})
    if not out:
        return None
    return out


def _compute_anomaly(metrics_before: Any, metrics_after: Any) -> bool:
    if not isinstance(metrics_before, dict) or not isinstance(metrics_after, dict):
        return False
    return any(metrics_before.get(key) != metrics_after.get(key) for key in _METRIC_NOOP_KEYS)


def _compute_noop(
    *,
    smiles_before: Any,
    smiles_after: Any,
    status: Any,
    error: Any,
    metrics_before: Any,
    metrics_after: Any,
    end_missing: bool,
) -> bool:
    if end_missing:
        return False
    if status != "ok" or error is not None:
        return False
    if not isinstance(smiles_before, str) or not isinstance(smiles_after, str):
        return False
    if smiles_before.strip() == "" or smiles_after.strip() == "":
        return False
    if smiles_before != smiles_after:
        return False
    return _same_metrics(metrics_before, metrics_after)


def _compute_always_visible(*, event_type: Any, status: Any, error: Any, meta: Any) -> bool:
    if event_type == "trace_truncated":
        return True
    if status != "ok":
        return True
    if error is not None:
        return True
    return bool(isinstance(meta, dict) and any(key in meta for key in _ALWAYS_VISIBLE_META_KEYS))


def _build_ui_nodes(
    events: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    span_end_candidates: dict[tuple[int, str], list[tuple[int, dict[str, Any]]]] = {}
    for idx, event in enumerate(events):
        if event.get("event_type") != "span_end":
            continue
        parent_id = event.get("parent_id")
        raw_op = event.get("op")
        if isinstance(parent_id, int) and isinstance(raw_op, str):
            span_end_candidates.setdefault((parent_id, raw_op), []).append((idx, event))

    nodes: list[dict[str, Any]] = []
    id_to_node: dict[int, dict[str, Any]] = {}
    used_end_event_ids: set[int] = set()

    for idx, event in enumerate(events):
        event_type = event.get("event_type")
        event_id = event.get("event_id")
        if not isinstance(event_id, int):
            continue
        if event_type == "span_end":
            continue

        raw_parent_id = event.get("parent_id")
        parent_id = raw_parent_id if isinstance(raw_parent_id, int) else None
        raw_op = event.get("op")
        op: str = raw_op if isinstance(raw_op, str) else ""
        raw_phase = event.get("phase")
        phase: str = raw_phase if isinstance(raw_phase, str) else ""
        status = event.get("status")
        error = event.get("error")
        smiles_before = event.get("smiles_in")
        smiles_after = event.get("smiles_out")
        metrics_before: Any = None
        metrics_after: Any = None
        end_event: Optional[dict[str, Any]] = None
        end_missing = False
        orphan_span_ends: list[dict[str, Any]] = []

        if event_type == "span_start":
            metrics_before = event.get("metrics")
            candidate_key = (event_id, op)
            candidates = span_end_candidates.get(candidate_key, [])
            for end_idx, candidate in candidates:
                candidate_id = candidate.get("event_id")
                if not isinstance(candidate_id, int):
                    continue
                if end_idx <= idx or candidate_id in used_end_event_ids:
                    continue
                end_event = candidate
                used_end_event_ids.add(candidate_id)
                break

            if end_event is None:
                end_missing = True
                status = event.get("status")
                error = event.get("error")
                smiles_before = event.get("smiles_in")
                smiles_after = None
                metrics_after = None
                meta = _merge_meta(event.get("meta"), None)
            else:
                status = end_event.get("status", event.get("status"))
                error = end_event.get("error", event.get("error"))
                smiles_before = event.get("smiles_in")
                smiles_after = end_event.get("smiles_out")
                metrics_after = end_event.get("metrics")
                meta = _merge_meta(event.get("meta"), end_event.get("meta"))

                end_idx = next(
                    (
                        candidate_idx
                        for candidate_idx, candidate in candidates
                        if candidate is end_event
                    ),
                    idx,
                )
                for extra_idx, extra in candidates:
                    extra_id = extra.get("event_id")
                    if not isinstance(extra_id, int):
                        continue
                    if extra_idx <= end_idx:
                        continue
                    if extra_id in used_end_event_ids:
                        continue
                    orphan_span_ends.append(extra)
        else:
            meta = _merge_meta(event.get("meta"), None)

        noop = _compute_noop(
            smiles_before=smiles_before,
            smiles_after=smiles_after,
            status=status,
            error=error,
            metrics_before=metrics_before,
            metrics_after=metrics_after,
            end_missing=end_missing,
        )
        always_visible = _compute_always_visible(
            event_type=event_type,
            status=status,
            error=error,
            meta=meta,
        )
        lowered_op = op.lower()
        if lowered_op.endswith("validate_omol") or lowered_op.endswith("omol_score"):
            noop = False
            always_visible = True

        smiles_changed = False
        if isinstance(meta, dict):
            smiles_changed = bool(meta.get("changed"))
        anomaly = _compute_anomaly(metrics_before, metrics_after)
        delta_table = _metric_delta_table(metrics_before, metrics_after)
        changed = (
            bool(smiles_changed)
            or bool((status != "ok") or (error is not None))
            or bool(end_missing)
        )

        node: dict[str, Any] = {
            "index": idx,
            "node_id": event_id,
            "parent_id": parent_id,
            "event_type": event_type,
            "op": op,
            "phase": phase,
            "status": status,
            "error": error,
            "smiles_before": smiles_before,
            "smiles_after": smiles_after,
            "meta": meta,
            "metrics_before": metrics_before,
            "metrics_after": metrics_after,
            "delta_table": delta_table,
            "noop": noop,
            "always_visible": always_visible,
            "is_success": status == "ok",
            "is_failure": (status != "ok") or (error is not None),
            "changed": changed,
            "anomaly": anomaly,
            "smiles_changed": smiles_changed,
            "is_final_used_path": False,
            "is_final_used_result": False,
            "end_missing": end_missing,
            "end_event": end_event,
            "orphan_span_ends": orphan_span_ends,
            "event": event,
            "children": [],
        }
        nodes.append(node)
        if event_id not in id_to_node:
            id_to_node[event_id] = node

    roots: list[dict[str, Any]] = []
    for node in nodes:
        parent_id = node.get("parent_id")
        if isinstance(parent_id, int) and parent_id in id_to_node:
            parent = id_to_node[parent_id]
            if parent is not node:
                parent["children"].append(node)
                continue
        roots.append(node)

    orphan_span_ends: list[dict[str, Any]] = []
    for _, event in enumerate(events):
        if event.get("event_type") != "span_end":
            continue
        event_id = event.get("event_id")
        if not isinstance(event_id, int):
            orphan_span_ends.append(event)
            continue
        if event_id not in used_end_event_ids:
            orphan_span_ends.append(event)

    return roots, orphan_span_ends


def _render_smiles_block(renderer: _SmilesRenderer, label: str, smiles: Optional[str]) -> str:
    text = html.escape(smiles) if isinstance(smiles, str) else ""
    svg = renderer.render(smiles)
    if svg is not None:
        body = f'<div class="mol-svg">{svg}</div><pre class="smiles-text">{text or "(none)"}</pre>'
    else:
        body = f'<pre class="smiles-text">{text or "(none)"}</pre>'
    return f'<section class="smiles"><h4>{label}</h4>{body}</section>'


def _node_summary(node: dict[str, Any]) -> tuple[str, str]:
    event_id_raw = str(node.get("node_id"))
    parent_id_raw = str(node.get("parent_id"))
    event_type_raw = str(node.get("event_type", ""))
    op_raw = str(node.get("op", ""))
    phase_raw = str(node.get("phase", ""))
    status_raw = str(node.get("status", ""))
    smiles_in = node.get("smiles_before")
    smiles_out = node.get("smiles_after")

    event_id = html.escape(event_id_raw)
    parent_id = html.escape(parent_id_raw)
    event_type = html.escape(event_type_raw)
    op = html.escape(op_raw)
    phase = html.escape(phase_raw)
    status = html.escape(status_raw)

    summary = f"#{event_id} {op} [{phase}/{status}] type={event_type} parent={parent_id}"
    search_text = " ".join(
        [
            event_id_raw,
            op_raw,
            phase_raw,
            status_raw,
            event_type_raw,
            parent_id_raw,
            smiles_in if isinstance(smiles_in, str) else "",
            smiles_out if isinstance(smiles_out, str) else "",
        ]
    ).lower()
    return summary, search_text


def _mark_visible_by_default(nodes: list[dict[str, Any]]) -> None:
    def walk(node: dict[str, Any]) -> bool:
        child_visible = False
        for child in node["children"]:
            child_visible = walk(child) or child_visible
        visible = bool(node.get("always_visible")) or (not bool(node.get("noop"))) or child_visible
        node["visible_by_default"] = visible
        return visible

    for root in nodes:
        walk(root)


def _mark_final_used_path(nodes: list[dict[str, Any]]) -> None:
    id_to_node: dict[int, dict[str, Any]] = {}
    for node in nodes:
        node_id = node.get("node_id")
        if isinstance(node_id, int) and node_id not in id_to_node:
            id_to_node[node_id] = node
        node["is_final_used_path"] = False
        node["is_final_used_result"] = False

    target: Optional[dict[str, Any]] = None
    for node in reversed(nodes):
        if node.get("smiles_after") is not None:
            target = node
            break

    if target is None:
        return

    target["is_final_used_result"] = True

    visited: set[int] = set()
    current: Optional[dict[str, Any]] = target
    while current is not None:
        current["is_final_used_path"] = True
        node_id = current.get("node_id")
        if not isinstance(node_id, int):
            break
        if node_id in visited:
            break
        visited.add(node_id)

        parent_id = current.get("parent_id")
        if not isinstance(parent_id, int):
            break
        current = id_to_node.get(parent_id)


def _mark_scoring_candidate_path(nodes: list[dict[str, Any]]) -> bool:
    id_to_node: dict[int, dict[str, Any]] = {}
    candidate_smiles: set[str] = set()

    for node in nodes:
        node_id = node.get("node_id")
        if isinstance(node_id, int) and node_id not in id_to_node:
            id_to_node[node_id] = node
        node["is_scoring_candidate_path"] = False

        meta = node.get("meta")
        if isinstance(meta, dict):
            raw_candidates = meta.get("candidate_smiles")
            if isinstance(raw_candidates, str):
                text = raw_candidates.strip()
                if text != "":
                    candidate_smiles.add(text)
            elif isinstance(raw_candidates, list):
                for item in raw_candidates:
                    if not isinstance(item, str):
                        continue
                    text = item.strip()
                    if text != "":
                        candidate_smiles.add(text)

        op = node.get("op")
        smiles_before = node.get("smiles_before")
        if (
            isinstance(op, str)
            and "omol_score" in op.lower()
            and isinstance(smiles_before, str)
            and smiles_before.strip() != ""
        ):
            candidate_smiles.add(smiles_before)

    if not candidate_smiles:
        return False

    candidate_nodes: list[dict[str, Any]] = []
    for node in nodes:
        smiles_before = node.get("smiles_before")
        smiles_after = node.get("smiles_after")
        if isinstance(smiles_before, str) and smiles_before in candidate_smiles:
            candidate_nodes.append(node)
            continue
        if isinstance(smiles_after, str) and smiles_after in candidate_smiles:
            candidate_nodes.append(node)

    for node in candidate_nodes:
        current: Optional[dict[str, Any]] = node
        visited: set[int] = set()
        while current is not None:
            current["is_scoring_candidate_path"] = True
            node_id = current.get("node_id")
            if not isinstance(node_id, int):
                break
            if node_id in visited:
                break
            visited.add(node_id)

            parent_id = current.get("parent_id")
            if not isinstance(parent_id, int):
                break
            current = id_to_node.get(parent_id)

    return True


def _build_graph_model(
    nodes: list[dict[str, Any]], *, scoring_candidate_active: bool
) -> dict[str, Any]:
    graph_nodes: list[dict[str, Any]] = []
    graph_edges: list[dict[str, Any]] = []

    node_ids = {
        node_id for node in nodes for node_id in [node.get("node_id")] if isinstance(node_id, int)
    }

    for node in nodes:
        node_id = node.get("node_id")
        parent_id = node.get("parent_id")
        graph_nodes.append(
            {
                "id": node_id,
                "parent_id": parent_id,
                "label": node.get("op"),
                "status": node.get("status"),
                "error": node.get("error"),
                "noop": bool(node.get("noop")),
                "always_visible": bool(node.get("always_visible")),
                "is_success": bool(node.get("is_success")),
                "is_failure": bool(node.get("is_failure")),
                "is_final_used_path": bool(node.get("is_final_used_path")),
                "is_final_used_result": bool(node.get("is_final_used_result")),
                "is_scoring_candidate_path": bool(node.get("is_scoring_candidate_path")),
            }
        )
        if isinstance(node_id, int) and isinstance(parent_id, int) and parent_id in node_ids:
            graph_edges.append({"source": parent_id, "target": node_id})

    return {
        "nodes": graph_nodes,
        "edges": graph_edges,
        "scoring_candidate_active": bool(scoring_candidate_active),
    }


def _build_candidate_pool(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidate_pool: list[dict[str, Any]] = []
    for node in nodes:
        op = node.get("op")
        if not isinstance(op, str) or "omol_score" not in op.lower():
            continue

        node_id = node.get("node_id")
        if not isinstance(node_id, int):
            continue

        smiles = node.get("smiles_before")
        if not isinstance(smiles, str):
            continue

        meta = node.get("meta")
        if not isinstance(meta, dict):
            continue
        score = meta.get("omol_score")
        if not isinstance(score, (int, float)):
            continue

        candidate_pool.append(
            {
                "node_id": node_id,
                "smiles": smiles,
                "score": float(score),
            }
        )

    candidate_pool.sort(key=lambda candidate: candidate["score"])
    return candidate_pool


def _node_data_attrs(node: dict[str, Any], *, search_text: str) -> str:
    attrs = [
        f'data-node-id="{html.escape(str(node.get("node_id")), quote=True)}"',
        f'data-filter-text="{html.escape(search_text, quote=True)}"',
    ]
    if node.get("noop"):
        attrs.append('data-noop="1"')
    if node.get("always_visible"):
        attrs.append('data-always-visible="1"')
    if node.get("is_success"):
        attrs.append('data-is-success="1"')
    if node.get("is_failure"):
        attrs.append('data-is-failure="1"')
    if node.get("is_final_used_path"):
        attrs.append('data-is-final-used-path="1"')
    if node.get("is_scoring_candidate_path"):
        attrs.append('data-is-scoring-candidate-path="1"')
    if node.get("changed"):
        attrs.append('data-changed="1"')
    if node.get("anomaly"):
        attrs.append('data-anomaly="1"')
    if not bool(node.get("visible_by_default", True)):
        attrs.append('data-default-hidden="1"')
    return " ".join(attrs)


def _encode_vscode_path(path: str, *, keep_windows_drive_colon: bool) -> str:
    normalized = path.replace("\\", "/")
    has_leading_slash = normalized.startswith("/")
    segments = [segment for segment in normalized.split("/") if segment != ""]

    encoded_segments: list[str] = []
    for idx, segment in enumerate(segments):
        if (
            keep_windows_drive_colon
            and idx == 0
            and len(segment) == 2
            and segment[1] == ":"
            and segment[0].isalpha()
        ):
            encoded_segments.append(f"{segment[0]}:")
            continue
        encoded_segments.append(quote(segment, safe="._-~"))

    encoded = "/".join(encoded_segments)
    if has_leading_slash:
        return f"/{encoded}" if encoded else "/"
    return encoded


def _build_vscode_source_links(src_path: str, line: int, col: int = 1) -> dict[str, str]:
    normalized_path = src_path.replace("\\", "/")
    line_num = line if line > 0 else 1
    col_num = col if col > 0 else 1

    encoded_local = _encode_vscode_path(normalized_path, keep_windows_drive_colon=True)
    fallback_href = f"vscode://file/{encoded_local}:{line_num}:{col_num}"

    distro = os.environ.get("WSL_DISTRO_NAME")
    primary_href = fallback_href
    fallback_command = f'code --goto "{normalized_path}:{line_num}:{col_num}"'
    if isinstance(distro, str) and distro.strip() != "":
        distro_name = distro.strip()
        encoded_distro = quote(distro_name, safe="._-~")
        encoded_wsl_path = _encode_vscode_path(normalized_path, keep_windows_drive_colon=False)
        if not encoded_wsl_path.startswith("/"):
            encoded_wsl_path = f"/{encoded_wsl_path}"
        primary_href = (
            f"vscode://vscode-remote/wsl+{encoded_distro}{encoded_wsl_path}:{line_num}:{col_num}"
        )
        fallback_command = (
            f'code --remote wsl+{distro_name} --goto "{normalized_path}:{line_num}:{col_num}"'
        )

    return {
        "primary_href": primary_href,
        "fallback_href": fallback_href,
        "fallback_command": fallback_command,
    }


def _render_source_link(meta: dict[str, Any]) -> str:
    src_path = meta.get("src_path")
    src_line = meta.get("src_line")
    src_qualname = meta.get("src_qualname")
    src_module = meta.get("src_module")
    if not isinstance(src_path, str) or src_path.strip() == "":
        return ""
    line = int(src_line) if isinstance(src_line, int) and src_line > 0 else 1
    display = f"{src_path}:{line}"
    if isinstance(src_qualname, str) and src_qualname.strip() != "":
        display = f"{display} ({src_qualname})"
    display_escaped = html.escape(display)
    module_note = ""
    if isinstance(src_module, str) and src_module.strip() != "":
        module_note = f"<span> {html.escape(src_module)}</span>"

    source_links = _build_vscode_source_links(src_path, line, 1)
    primary_href_escaped = html.escape(source_links["primary_href"], quote=True)
    fallback_href_escaped = html.escape(source_links["fallback_href"], quote=True)
    fallback_command_escaped = html.escape(source_links["fallback_command"])
    open_label = (
        "open (wsl)" if source_links["primary_href"] != source_links["fallback_href"] else "open"
    )

    return " ".join(
        [
            f'<a href="{primary_href_escaped}">{display_escaped}</a>{module_note}',
            f'<a href="{primary_href_escaped}">{html.escape(open_label)}</a>',
            f'<a href="{fallback_href_escaped}">open (file)</a>',
            f"<code>{fallback_command_escaped}</code>",
        ]
    )


def _render_delta_table(delta_table: Any) -> str:
    if not isinstance(delta_table, list) or not delta_table:
        return ""
    rows: list[str] = []
    for row in delta_table:
        if not isinstance(row, dict):
            continue
        metric = html.escape(str(row.get("metric", "")))
        before = html.escape(str(row.get("before", "")))
        after = html.escape(str(row.get("after", "")))
        delta = row.get("delta")
        delta_text = ""
        delta_text = f"{delta:+g}" if isinstance(delta, (int, float)) else str(row.get("delta", ""))
        rows.append(
            f"<tr><td>{metric}</td><td>{before}</td><td>{after}</td><td>{html.escape(delta_text)}</td></tr>"
        )
    if not rows:
        return ""
    return "".join(
        [
            '<section class="delta-table">',
            "<h4>metrics delta</h4>",
            '<table class="delta-table-grid">',
            "<thead><tr><th>metric</th><th>before</th><th>after</th><th>delta</th></tr></thead>",
            f"<tbody>{''.join(rows)}</tbody>",
            "</table>",
            "</section>",
        ]
    )


def _render_tree_node(node: dict[str, Any], *, depth: int, seen: set[int]) -> str:
    node_idx = int(node["index"])
    if node_idx in seen:
        return '<li class="cycle">cycle detected</li>'
    seen.add(node_idx)

    node_id = html.escape(str(node.get("node_id")))
    summary, search_text = _node_summary(node)
    data_attrs = _node_data_attrs(node, search_text=search_text)
    hidden_attr = " hidden" if not bool(node.get("visible_by_default", True)) else ""
    link = (
        '<a class="tree-link" '
        f'data-depth="{depth}" '
        f"{data_attrs}{hidden_attr} "
        f'href="#event-{node_id}">{html.escape(summary)}</a>'
    )

    parts = [f'<li class="tree-item" data-node-id="{node_id}">', link]
    children = node["children"]
    if children:
        parts.append('<div class="children"><ul class="tree-list">')
        for child in children:
            parts.append(_render_tree_node(child, depth=depth + 1, seen=seen))
        parts.append("</ul></div>")
    parts.append("</li>")

    seen.remove(node_idx)
    return "".join(parts)


def _render_detail_card(node: dict[str, Any], *, renderer: _SmilesRenderer) -> str:
    event_id_raw = str(node.get("node_id"))
    parent_id_raw = str(node.get("parent_id"))
    event_type_raw = str(node.get("event_type", ""))
    op_raw = str(node.get("op", ""))
    phase_raw = str(node.get("phase", ""))
    status_raw = str(node.get("status", ""))
    smiles_in = node.get("smiles_before")
    smiles_out = node.get("smiles_after")
    noop_raw = str(bool(node.get("noop"))).lower()
    always_visible_raw = str(bool(node.get("always_visible"))).lower()
    end_missing_raw = str(bool(node.get("end_missing"))).lower()

    event_id = html.escape(event_id_raw)
    parent_id = html.escape(parent_id_raw)
    event_type = html.escape(event_type_raw)
    op = html.escape(op_raw)
    phase = html.escape(phase_raw)
    status = html.escape(status_raw)
    summary, search_text = _node_summary(node)
    data_attrs = _node_data_attrs(node, search_text=search_text)
    hidden_attr = " hidden" if not bool(node.get("visible_by_default", True)) else ""

    meta = node.get("meta")
    source_html = ""
    if isinstance(meta, dict):
        source_html = _render_source_link(meta)

    delta_html = _render_delta_table(node.get("delta_table"))

    return "".join(
        [
            f'<article class="detail" id="event-{event_id}" {data_attrs}{hidden_attr}>',
            f"<h3>{html.escape(summary)}</h3>",
            '<dl class="fields">',
            f"<dt>event_id</dt><dd>{event_id}</dd>",
            f"<dt>parent_id</dt><dd>{parent_id}</dd>",
            f"<dt>event_type</dt><dd>{event_type}</dd>",
            f"<dt>op</dt><dd>{op}</dd>",
            f"<dt>phase</dt><dd>{phase}</dd>",
            f"<dt>status</dt><dd>{status}</dd>",
            f"<dt>noop</dt><dd>{html.escape(noop_raw)}</dd>",
            f"<dt>always_visible</dt><dd>{html.escape(always_visible_raw)}</dd>",
            f"<dt>end_missing</dt><dd>{html.escape(end_missing_raw)}</dd>",
            f"<dt>changed</dt><dd>{html.escape(str(bool(node.get('changed'))).lower())}</dd>",
            f"<dt>anomaly</dt><dd>{html.escape(str(bool(node.get('anomaly'))).lower())}</dd>",
            f"<dt>source</dt><dd>{source_html or ''}</dd>",
            "</dl>",
            delta_html,
            '<div class="smiles-row">',
            _render_smiles_block(renderer, "smiles_before", smiles_in),
            _render_smiles_block(renderer, "smiles_after", smiles_out),
            "</div>",
            "</article>",
        ]
    )


def _flatten_nodes(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def walk(current: dict[str, Any]) -> None:
        out.append(current)
        for child in current["children"]:
            walk(child)

    for node in nodes:
        walk(node)
    return out


def _template_env() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(
            enabled_extensions=("html", "htm", "xml", "j2"),
            default_for_string=True,
            default=True,
        ),
    )


def render_trace_report(
    trace_dir: Path,
    out_path: Optional[Path] = None,
    max_events_rendered: int = 2000,
    max_atoms_svg: int = 200,
) -> dict[str, Any]:
    trace_path = trace_dir / TRACE_FILENAME
    if out_path is None:
        out_path = trace_dir / REPORT_FILENAME

    events, valid_events, malformed_lines, truncated = _read_trace_events(
        trace_path,
        max_events_rendered=max_events_rendered,
    )

    roots, orphan_span_ends = _build_ui_nodes(events)
    _mark_visible_by_default(roots)
    renderer = _SmilesRenderer(max_atoms_svg=max_atoms_svg)
    tree_html: list[str] = []
    detail_html: list[str] = []
    seen: set[int] = set()

    flat_nodes = _flatten_nodes(roots)
    ordered_nodes = sorted(flat_nodes, key=lambda node: int(node["index"]))
    _mark_final_used_path(ordered_nodes)
    scoring_candidate_active = _mark_scoring_candidate_path(ordered_nodes)
    graph_model = _build_graph_model(
        ordered_nodes,
        scoring_candidate_active=scoring_candidate_active,
    )
    candidate_pool = _build_candidate_pool(ordered_nodes)
    graph_model_json = json.dumps(graph_model, separators=(",", ":")).replace("</", "<\\/")
    candidate_pool_json = json.dumps(candidate_pool, separators=(",", ":")).replace("</", "<\\/")
    for root in roots:
        tree_html.append(_render_tree_node(root, depth=0, seen=seen))
    for node in ordered_nodes:
        detail_html.append(_render_detail_card(node, renderer=renderer))

    note = ""
    if truncated:
        note = (
            '<p class="note">Truncated view: showing first '
            f"{len(events)} events out of {valid_events}.</p>"
        )

    orphan_note = ""
    if orphan_span_ends:
        orphan_note = (
            '<p class="note">Orphan span_end events retained for later rendering: '
            f"{len(orphan_span_ends)}.</p>"
        )

    template = _template_env().get_template(_TEMPLATE_NAME)
    html_text = template.render(
        rendered_events=len(events),
        valid_events=valid_events,
        malformed_lines=malformed_lines,
        note_html=note,
        orphan_note_html=orphan_note,
        tree_html="".join(tree_html),
        detail_html="".join(detail_html),
        graph_model_json=graph_model_json,
        candidate_pool_json=candidate_pool_json,
    )
    out_path.write_text(html_text, encoding="utf-8")

    return {
        "out_path": str(out_path),
        "trace_path": str(trace_path),
        "rendered_events": len(events),
        "valid_events": valid_events,
        "malformed_lines": malformed_lines,
        "truncated": truncated,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate HTML debug report from JSONL trace")
    parser.add_argument(
        "--trace-dir", type=Path, required=True, help="Directory containing trace.jsonl"
    )
    parser.add_argument("--out", type=Path, default=None, help="Output HTML path")
    parser.add_argument("--max-events-rendered", type=int, default=2000)
    parser.add_argument("--max-atoms-svg", type=int, default=200)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    result = render_trace_report(
        trace_dir=args.trace_dir,
        out_path=args.out,
        max_events_rendered=max(args.max_events_rendered, 1),
        max_atoms_svg=max(args.max_atoms_svg, 1),
    )
    print(result["out_path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
