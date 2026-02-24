from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional


class TraceWriter:
    def __init__(self, trace_path: Path) -> None:
        self._fh = trace_path.open("w", encoding="utf-8")
        self._next_event_id = 1
        self._stack: list[int] = []

    def close(self) -> None:
        self._fh.close()

    def _emit(self, event: dict[str, Any]) -> None:
        self._fh.write(json.dumps(event, ensure_ascii=True, separators=(",", ":")))
        self._fh.write("\n")

    def _new_event_id(self) -> int:
        event_id = self._next_event_id
        self._next_event_id += 1
        return event_id

    def _parent_id(self) -> Optional[int]:
        if self._stack:
            return self._stack[-1]
        return None

    def span_start(
        self,
        *,
        op: str,
        phase: str,
        smiles_in: Optional[str],
        meta: dict[str, Any],
        metrics: Optional[dict[str, Any]],
    ) -> int:
        event_id = self._new_event_id()
        event = {
            "event_type": "span_start",
            "event_id": event_id,
            "parent_id": self._parent_id(),
            "op": op,
            "phase": phase,
            "status": "start",
            "error": None,
            "smiles_in": smiles_in,
            "meta": meta,
            "metrics": metrics,
        }
        self._emit(event)
        self._stack.append(event_id)
        return event_id

    def span_end(
        self,
        *,
        start_event_id: int,
        op: str,
        status: str,
        error: Optional[str],
        smiles_out: Optional[str],
        meta: dict[str, Any],
        metrics: Optional[dict[str, Any]],
    ) -> None:
        if self._stack and self._stack[-1] == start_event_id:
            self._stack.pop()
        event = {
            "event_type": "span_end",
            "event_id": self._new_event_id(),
            "parent_id": start_event_id,
            "op": op,
            "status": status,
            "error": error,
            "smiles_out": smiles_out,
            "meta": meta,
            "metrics": metrics,
        }
        self._emit(event)
