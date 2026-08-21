from __future__ import annotations

import subprocess
from html.parser import HTMLParser
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]


class _IdTreeParser(HTMLParser):
    VOID_TAGS = {
        "area",
        "base",
        "br",
        "col",
        "embed",
        "hr",
        "img",
        "input",
        "link",
        "meta",
        "param",
        "source",
        "track",
        "wbr",
    }

    def __init__(self) -> None:
        super().__init__()
        self.stack: list[str | None] = []
        self.nodes: dict[str, tuple[str | None, list[str]]] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        node_id = values.get("id")
        parent_id = next((value for value in reversed(self.stack) if value), None)
        if node_id:
            self.nodes[node_id] = (parent_id, (values.get("class") or "").split())
        if tag not in self.VOID_TAGS:
            self.stack.append(node_id)

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        node_id = values.get("id")
        parent_id = next((value for value in reversed(self.stack) if value), None)
        if node_id:
            self.nodes[node_id] = (parent_id, (values.get("class") or "").split())

    def handle_endtag(self, tag: str) -> None:
        self.stack.pop()


def test_primary_structure_panels_are_fixed_and_legacy_switch_is_absent() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    parser = _IdTreeParser()
    parser.feed(html)

    expected = ["viewerPanel", "referenceXyzPanel", "primaryVisual", "secondaryVisual"]
    assert [node for node in expected if parser.nodes[node][0] == "primaryVisuals"] == expected
    assert "xyz-card" in parser.nodes["viewerPanel"][1]
    assert "reference-xyz-card" in parser.nodes["referenceXyzPanel"][1]
    assert "candidate-card" in parser.nodes["primaryVisual"][1]
    assert "reference-card" in parser.nodes["secondaryVisual"][1]
    assert html.index('id="viewerPanel"') < html.index('id="referenceXyzPanel"')
    assert html.index('id="referenceXyzPanel"') < html.index('id="primaryVisual"')
    assert html.index('id="primaryVisual"') < html.index('id="secondaryVisual"')
    assert "render-kind" not in html
    assert "primaryKind" not in javascript
    assert "secondaryKind" not in javascript
    assert 'querySelectorAll(".render-kind")' not in javascript
    assert 'loadRender("candidate", "primary", token)' in javascript
    assert 'loadRender("reference", "secondary", token)' in javascript


def test_review_history_ui_keeps_auto_next_and_restores_undone_case() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")

    assert 'id="undoReview"' in html
    assert 'id="lastReviewSummary"' in html
    assert 'id="recentReviewList"' in html
    assert "state.reviewHistory = state.reviewHistory.slice(0, 20)" in javascript
    assert "await loadCase(queuedNextCaseId)" in javascript
    assert "await loadCase(result.case_id)" in javascript
    assert '(event.ctrlKey || event.metaKey) && event.key.toLowerCase() === "z"' in javascript
    assert '"X-MolGR-Review-Session": reviewSessionId' in javascript


def test_optional_family_filter_hides_and_sidebar_counts_relocalize() -> None:
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    stylesheet = (APP_DIR / "static" / "style.css").read_text(encoding="utf-8")

    assert ".filters label[hidden] { display: none; }" in stylesheet
    toggle_language = javascript.split("function toggleLanguage()", 1)[1].split("\n}", 1)[0]
    assert "applyLanguage();" in toggle_language
    assert "loadStats().catch" in toggle_language


def test_mapped_comparison_has_no_persistent_atom_overlays() -> None:
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")

    assert "function applyMappedComparison" not in javascript
    assert "applyMappedComparison(viewer" not in javascript
    assert "function renderMappedComparisonNote" in javascript


def test_representative_reference_xyz_warns_without_hover_ui() -> None:
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")

    assert "data.mapping_is_representative" in javascript
    assert 'tr("representativeMappingWarning")' in javascript
    assert 'tr("ambiguityType")' in javascript
    assert 'tr("diagnosticReasonLabel")' in javascript
    assert "mappingRepresentativeHover" not in javascript


def test_reference_states_render_failures_and_stale_response_guard() -> None:
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    javascript = javascript.rsplit("\ninit().catch", 1)[0]
    prelude = r"""
globalThis.localStorage = { getItem() { return null; }, setItem() {} };
"""
    harness = r"""
const assert = (condition, message) => { if (!condition) throw new Error(message); };
const elements = new Map();
function element(id) {
  if (!elements.has(id)) {
    elements.set(id, {
      id, innerHTML: "", textContent: "", hidden: false, disabled: false, title: "",
      dataset: {}, classList: { toggle() {}, add() {}, remove() {} },
      querySelector(selector) { return selector.includes("svg") && this.innerHTML.includes("<svg") ? {} : null; },
      setAttribute() {}, removeAttribute() {},
    });
  }
  return elements.get(id);
}
globalThis.document = {
  getElementById: element,
  querySelectorAll() { return []; },
  documentElement: { style: { setProperty() {} } },
  body: element("body"),
};
globalThis.window = {
  location: { href: "http://localhost/" },
  history: { replaceState() {} },
  requestAnimationFrame(callback) { callback(); return 1; },
};

state.current = { case_id: "A", reference_smiles: "C", reference_parse_status: "ok" };
state.caseRequestToken = 1;
state.referenceRenderStatus = "render_failed";
assert(referenceState({ reference_smiles: "" }) === "missing", "missing state");
assert(referenceState({ reference_smiles: "bad", reference_parse_status: "error" }) === "parse_invalid", "parse-invalid state");
assert(referenceState(state.current) === "render_failed", "render-failed state");
state.referenceRenderStatus = "available";
assert(referenceState(state.current) === "available", "available state");

async function response(payload, status = 200) {
  return {
    ok: status >= 200 && status < 300,
    status,
    statusText: `HTTP ${status}`,
    headers: { get() { return "application/json"; } },
    async json() { return payload; },
  };
}

async function checkFailure(fetchImpl, expectedDetail) {
  state.current = { case_id: "A", reference_smiles: "C", reference_parse_status: "ok" };
  state.caseRequestToken += 1;
  state.referenceRenderStatus = "render_failed";
  state.referenceRenderError = "";
  element("primaryVisual").hidden = false;
  element("secondaryVisual").hidden = false;
  globalThis.fetch = fetchImpl;
  await loadRender("reference", "secondary", state.caseRequestToken);
  assert(referenceState(state.current) === "render_failed", "failure must remain render_failed");
  assert(!element("primaryVisual").hidden, "candidate panel hidden after reference failure");
  assert(!element("secondaryVisual").hidden, "reference panel hidden after reference failure");
  assert(element("secondarySvg").innerHTML.includes("Reference render failed"), "reviewer failure message");
  assert(state.referenceRenderError.includes(expectedDetail), `missing technical detail: ${expectedDetail}`);
}

(async () => {
  await checkFailure(async () => response({ error: "renderer rejected molecule" }), "payload.error=renderer rejected molecule");
  await checkFailure(async () => response({ error: "not found" }, 404), "HTTP 404");
  await checkFailure(async () => { throw new Error("network down"); }, "network down");
  await checkFailure(async () => response({ svg: "", smiles: "C" }), "empty SVG");

  state.current = { case_id: "A", reference_smiles: "C", reference_parse_status: "ok" };
  state.caseRequestToken += 1;
  globalThis.fetch = async () => response({ svg: "<svg></svg>", smiles: "C" });
  await loadRender("reference", "secondary", state.caseRequestToken);
  assert(referenceState(state.current) === "available", "successful SVG must be available");

  let releaseA;
  const delayed = new Promise((resolve) => { releaseA = resolve; });
  state.current = { case_id: "A", reference_smiles: "C", reference_parse_status: "ok" };
  state.caseRequestToken += 1;
  const tokenA = state.caseRequestToken;
  globalThis.fetch = async () => delayed;
  const pendingA = loadRender("reference", "secondary", tokenA);
  state.current = { case_id: "B", reference_smiles: "N", reference_parse_status: "ok" };
  state.caseRequestToken += 1;
  state.referenceRenderStatus = "render_failed";
  state.referenceRenderError = "B owns this state";
  element("secondarySvg").innerHTML = "B DOM";
  releaseA(await response({ svg: "<svg id='A'></svg>", smiles: "C" }));
  await pendingA;
  assert(state.current.case_id === "B", "stale request changed current case");
  assert(state.referenceRenderError === "B owns this state", "stale request changed reference state");
  assert(element("secondarySvg").innerHTML === "B DOM", "stale request changed B DOM");
})().catch((error) => { console.error(error.stack || error); process.exit(1); });
"""
    subprocess.run(
        ["node", "-e", f"{prelude}\n{javascript}\n{harness}"],
        check=True,
        cwd=APP_DIR,
        text=True,
    )
