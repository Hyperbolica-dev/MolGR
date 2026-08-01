const state = {
  cases: [],
  total: 0,
  offset: 0,
  limit: 80,
  current: null,
  viewStyle: "stick",
  primaryKind: "candidate",
  secondaryKind: "reference",
  xyzViewer: null,
  candidateViewer: null,
  currentCandidateSdf: "",
  currentLiveCandidate: null,
  ketcherLoaded: false,
  layout: null,
};

const labels = {
  candidate: "当前候选重建",
  reference: "Reference",
  candidate_organic: "候选图 organic",
  reference_organic: "Reference organic",
};

const statusLabels = {
  accept_candidate: "候选图为准",
  accept_reference: "参考图为准",
  accept_both: "接受两者",
  manual_reference: "人工修图为准",
  needs_followup: "待复核",
  skip: "跳过",
};

const layoutStorageKey = "moleculeReviewLayout.v1";
const layoutDefaults = {
  sidebarWidth: 320,
  viewerWidth: 420,
  reviewWidth: 440,
  ketcherHeight: 560,
};
const layoutMinimums = {
  sidebarWidth: 260,
  viewerWidth: 320,
  compareWidth: 420,
  reviewWidth: 360,
  ketcherHeight: 360,
};
const splitterSize = 10;
let viewerResizeFrame = 0;
let viewerResizeObserver = null;

function $(id) {
  return document.getElementById(id);
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(value, max));
}

function readStoredLayout() {
  try {
    const raw = localStorage.getItem(layoutStorageKey);
    if (!raw) return {};
    const parsed = JSON.parse(raw);
    return parsed && typeof parsed === "object" ? parsed : {};
  } catch (error) {
    return {};
  }
}

function saveLayout() {
  if (!state.layout) return;
  localStorage.setItem(layoutStorageKey, JSON.stringify(state.layout));
}

function queueViewerResize() {
  if (viewerResizeFrame) return;
  viewerResizeFrame = window.requestAnimationFrame(() => {
    viewerResizeFrame = 0;
    [state.xyzViewer, state.candidateViewer].forEach((viewer) => {
      if (viewer && typeof viewer.resize === "function") {
        viewer.resize();
        viewer.render();
      }
    });
  });
}

function applyLayout() {
  if (!state.layout) return;
  const root = document.documentElement.style;
  root.setProperty("--sidebar-width", `${state.layout.sidebarWidth}px`);
  root.setProperty("--viewer-width", `${state.layout.viewerWidth}px`);
  root.setProperty("--review-width", `${state.layout.reviewWidth}px`);
  root.setProperty("--ketcher-height", `${state.layout.ketcherHeight}px`);
  queueViewerResize();
}

function mainLayoutWidth() {
  const layout = $("mainLayout");
  if (!layout) return 0;
  return layout.getBoundingClientRect().width;
}

function clampLayoutValue(key, proposed, snapshot = state.layout || layoutDefaults) {
  if (key === "sidebarWidth") {
    const max = Math.max(layoutMinimums.sidebarWidth, window.innerWidth - 720);
    return clamp(proposed, layoutMinimums.sidebarWidth, max);
  }
  if (key === "viewerWidth") {
    const total = mainLayoutWidth();
    const max = Math.max(
      layoutMinimums.viewerWidth,
      total - layoutMinimums.compareWidth - splitterSize,
    );
    return clamp(proposed, layoutMinimums.viewerWidth, max);
  }
  if (key === "reviewWidth") {
    const total = mainLayoutWidth();
    const max = Math.max(
      layoutMinimums.reviewWidth,
      total - snapshot.viewerWidth - layoutMinimums.compareWidth - splitterSize * 2,
    );
    return clamp(proposed, layoutMinimums.reviewWidth, max);
  }
  if (key === "ketcherHeight") {
    const max = Math.max(layoutMinimums.ketcherHeight, window.innerHeight * 2);
    return clamp(proposed, layoutMinimums.ketcherHeight, max);
  }
  return proposed;
}

function initializeLayout() {
  state.layout = { ...layoutDefaults, ...readStoredLayout() };
  state.layout.sidebarWidth = clampLayoutValue("sidebarWidth", state.layout.sidebarWidth);
  state.layout.viewerWidth = clampLayoutValue("viewerWidth", state.layout.viewerWidth, state.layout);
  state.layout.reviewWidth = clampLayoutValue("reviewWidth", state.layout.reviewWidth, state.layout);
  state.layout.ketcherHeight = clampLayoutValue("ketcherHeight", state.layout.ketcherHeight);
  applyLayout();
}

function applyLayoutDelta(key, delta, snapshot) {
  if (!state.layout) return;
  let nextValue = snapshot[key];
  if (key === "reviewWidth") {
    nextValue = snapshot.reviewWidth - delta;
  } else {
    nextValue = snapshot[key] + delta;
  }
  state.layout[key] = clampLayoutValue(key, nextValue, snapshot);
  applyLayout();
}

function updateBodyResizeCursor(key, active) {
  document.body.classList.toggle("resizing", active);
  document.body.classList.toggle("resizing-row", active && key === "ketcherHeight");
}

function setupResizer(handle, key, axis) {
  if (!handle) return;
  const step = key === "ketcherHeight" ? 40 : 24;
  handle.addEventListener("pointerdown", (event) => {
    if (window.matchMedia("(max-width: 1120px)").matches) return;
    if (!state.layout) return;
    event.preventDefault();
    const snapshot = { ...state.layout };
    const origin = axis === "x" ? event.clientX : event.clientY;
    handle.classList.add("dragging");
    updateBodyResizeCursor(key, true);

    const onMove = (moveEvent) => {
      const delta = (axis === "x" ? moveEvent.clientX : moveEvent.clientY) - origin;
      applyLayoutDelta(key, delta, snapshot);
    };
    const stop = () => {
      window.removeEventListener("pointermove", onMove);
      window.removeEventListener("pointerup", stop);
      window.removeEventListener("pointercancel", stop);
      handle.classList.remove("dragging");
      updateBodyResizeCursor(key, false);
      saveLayout();
    };

    window.addEventListener("pointermove", onMove);
    window.addEventListener("pointerup", stop);
    window.addEventListener("pointercancel", stop);
  });

  handle.addEventListener("keydown", (event) => {
    if (!state.layout) return;
    const keyName = event.key;
    const isNegative = (axis === "x" && keyName === "ArrowLeft") || (axis === "y" && keyName === "ArrowUp");
    const isPositive = (axis === "x" && keyName === "ArrowRight") || (axis === "y" && keyName === "ArrowDown");
    if (!isNegative && !isPositive) return;
    event.preventDefault();
    const snapshot = { ...state.layout };
    const delta = (isNegative ? -1 : 1) * step;
    applyLayoutDelta(key, delta, snapshot);
    saveLayout();
  });
}

function bindLayoutResizers() {
  setupResizer(document.querySelector('[data-resize="sidebar"]'), "sidebarWidth", "x");
  setupResizer(document.querySelector('[data-resize="ketcher"]'), "ketcherHeight", "y");
  window.addEventListener(
    "resize",
    debounce(() => {
      if (!state.layout) return;
      state.layout.sidebarWidth = clampLayoutValue("sidebarWidth", state.layout.sidebarWidth);
      state.layout.viewerWidth = clampLayoutValue("viewerWidth", state.layout.viewerWidth, state.layout);
      state.layout.reviewWidth = clampLayoutValue("reviewWidth", state.layout.reviewWidth, state.layout);
      state.layout.ketcherHeight = clampLayoutValue("ketcherHeight", state.layout.ketcherHeight);
      applyLayout();
      saveLayout();
    }, 80),
  );
}

function observeViewerContainer() {
  if (!window.ResizeObserver) return;
  const containers = [$("viewer3d"), $("viewerCandidate3d")].filter(Boolean);
  if (!containers.length) return;
  viewerResizeObserver = new ResizeObserver(() => queueViewerResize());
  containers.forEach((container) => viewerResizeObserver.observe(container));
}

async function api(path, options = {}) {
  const response = await fetch(path, {
    headers: { "Content-Type": "application/json", ...(options.headers || {}) },
    ...options,
  });
  const contentType = response.headers.get("content-type") || "";
  const payload = contentType.includes("application/json") ? await response.json() : await response.text();
  if (!response.ok) {
    const message = typeof payload === "string" ? payload : payload.error || response.statusText;
    throw new Error(message);
  }
  return payload;
}

function badge(text, kind = "") {
  return `<span class="badge ${kind}">${escapeHtml(text || "")}</span>`;
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function setImageZoomState(box, label) {
  const image = box.querySelector("svg, img");
  box.classList.toggle("is-zoomable", Boolean(image));
  if (!image) {
    box.removeAttribute("role");
    box.removeAttribute("tabindex");
    box.removeAttribute("aria-label");
    delete box.dataset.zoomLabel;
    return;
  }
  box.setAttribute("role", "button");
  box.tabIndex = 0;
  box.setAttribute("aria-label", `放大${label}`);
  box.dataset.zoomLabel = label;
}

function openImageLightbox(box) {
  const image = box.querySelector("svg, img");
  if (!image) return;
  const dialog = $("imageLightbox");
  $("imageLightboxTitle").textContent = box.dataset.zoomLabel || "结构图";
  $("imageLightboxContent").replaceChildren(image.cloneNode(true));
  if (!dialog.open) dialog.showModal();
}

function closeImageLightbox() {
  const dialog = $("imageLightbox");
  if (dialog.open) dialog.close();
}

function categoryKind(category) {
  if (category === "graph_not_equivalent") return "warn";
  if (category === "backend_mismatch") return "bad";
  if (category === "python_version_mismatch") return "warn";
  if (category === "reference_not_comparable") return "warn";
  if (category === "reference_formula_mismatch") return "bad";
  if (category === "candidate_failed") return "bad";
  return "";
}

async function loadStats() {
  const stats = await api("/api/stats");
  const categoryLabels = {
    graph_not_equivalent: "图不等价",
    missing_reference_smiles: "缺参考",
    candidate_failed: "候选生成失败",
    backend_mismatch: "后端分歧",
    python_version_mismatch: "版本分歧",
    reference_not_comparable: "参考不可比",
    reference_formula_mismatch: "参考氢数不守恒",
  };
  const categoryEntries = Object.entries(stats.categories || {})
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([key, value]) => [categoryLabels[key] || key, value]);
  const entries = [
    ["全部", Object.values(stats.categories || {}).reduce((a, b) => a + b, 0)],
    ["未审核", (stats.review_statuses || {}).unreviewed || 0],
    ["候选图为准", (stats.review_statuses || {}).accept_candidate || 0],
    ["参考图为准", (stats.review_statuses || {}).accept_reference || 0],
    ["接受两者", (stats.review_statuses || {}).accept_both || 0],
    ["人工修图", (stats.review_statuses || {}).manual_reference || 0],
    ...categoryEntries,
  ];
  $("stats").innerHTML = entries
    .map(([name, value]) => `<div class="stat"><b>${value}</b>${escapeHtml(name)}</div>`)
    .join("");
  const metadata = stats.metadata || {};
  const runtime = stats.runtime || {};
  const metaEntries = [
    ["source", metadata.source_csv || ""],
    ["cases", metadata.record_count || ""],
    ["imported_at", metadata.imported_at || ""],
    ["checkout", `${runtime.git_revision || "unknown"}${runtime.git_dirty ? " (dirty)" : ""}`],
    ["python", runtime.python || ""],
    ["cpp", runtime.cpp_extension || ""],
  ].filter(([, value]) => value);
  $("datasetMeta").innerHTML = metaEntries
    .map(
      ([name, value]) =>
        `<div><span>${escapeHtml(name)}</span>` +
        `<code title="${escapeHtml(value)}">${escapeHtml(value)}</code></div>`,
    )
    .join("");
}

async function loadCases(reset = false) {
  if (reset) state.offset = 0;
  const params = new URLSearchParams();
  const category = $("categoryFilter").value;
  const status = $("statusFilter").value;
  const q = $("searchBox").value.trim();
  if (category) params.set("category", category);
  if (status) params.set("status", status);
  if (q) params.set("q", q);
  params.set("limit", state.limit);
  params.set("offset", state.offset);
  const data = await api(`/api/cases?${params.toString()}`);
  state.cases = data.items || [];
  state.total = data.total || 0;
  renderCaseList();
}

function renderCaseList() {
  $("caseList").innerHTML = state.cases
    .map((item) => {
      const selected = state.current && state.current.case_id === item.case_id ? "selected" : "";
      const status = item.review_status ? statusLabels[item.review_status] || item.review_status : "未审核";
      const fixture = item.fixture ? `fixture · ${item.fixture.kind}` : "无 fixture";
      return `
        <button class="case-item ${selected}" data-case-id="${escapeHtml(item.case_id)}" type="button">
          <span class="row"><strong>${escapeHtml(item.case_id)}</strong><span>#${item.row_index}</span></span>
          <span class="row">
            ${badge(item.category, categoryKind(item.category))}
            ${badge(status, item.review_status ? "ok" : "")}
            ${badge(fixture, item.fixture ? "ok" : "warn")}
          </span>
        </button>`;
    })
    .join("");
  $("pageInfo").textContent = `${state.offset + 1}-${Math.min(state.offset + state.limit, state.total)} / ${state.total}`;
  document.querySelectorAll(".case-item").forEach((button) => {
    button.addEventListener("click", () => loadCase(button.dataset.caseId));
  });
}

async function loadCase(caseId) {
  const item = await api(`/api/cases/${encodeURIComponent(caseId)}`);
  state.current = item;
  const url = new URL(window.location.href);
  url.searchParams.set("case", item.case_id);
  window.history.replaceState({}, "", url);
  state.primaryKind = "candidate";
  state.secondaryKind = "reference";
  state.currentLiveCandidate = null;
  renderCaseHeader();
  populateReviewForm();
  renderVersionComparison();
  renderDiagnostics();
  await Promise.all([loadXyz(item), loadCandidateSdf(item)]);
  await loadPair();
  renderCaseList();
}

function renderCaseHeader() {
  const item = state.current;
  const openTrace = $("openTrace");
  openTrace.disabled = !item;
  openTrace.title = item ? `在新窗口打开 ${item.case_id} 的 Trace` : "请先选择 case";
  const fixture = item.fixture;
  const fixtureBadge = fixture
    ? badge(`fixture · ${fixture.kind}`, "ok")
    : badge("无 fixture", "warn");
  $("caseTitle").textContent = `${item.case_id} · row ${item.row_index}`;
  const liveMismatchBadge = item.live_matches_candidate_snapshot === false
    ? badge("当前重建 != 候选快照", "bad")
    : "";
  $("caseMeta").innerHTML = [
    badge(item.category, categoryKind(item.category)),
    badge(
      `候选快照 · ${item.candidate_snapshot_status || "missing"}`,
      item.candidate_snapshot_status === "ok" ? "ok" : "bad",
    ),
    item.review_status ? badge(statusLabels[item.review_status] || item.review_status, "ok") : badge("未审核"),
    fixtureBadge,
    liveMismatchBadge,
  ].join(" ");
  const removeFixture = $("removeFixture");
  removeFixture.hidden = !fixture;
  removeFixture.disabled = !fixture;
  removeFixture.title = fixture
    ? `移除 ${fixture.structure_file} 并将审核状态设为待复核`
    : "当前 case 没有 fixture";
}

function populateReviewForm() {
  const item = state.current;
  $("correctedSmiles").value = item.corrected_smiles || "";
  $("correctedMolblock").value = item.corrected_molblock || "";
  $("notes").value = item.notes || "";
  $("reviewer").value = item.reviewer || localStorage.getItem("moleculeReviewReviewer") || "";
  document.querySelectorAll(".decision").forEach((button) => {
    button.classList.toggle("selected", button.dataset.status === item.review_status);
  });
}

function renderDiagnostics() {
  const item = state.current;
  const pairs = [
    ["total_charge", item.total_charge],
    ["spin_multiplicity", item.spin_multiplicity],
    ["total_radical_electrons", item.total_radical_electrons],
    ["reference_smiles", item.reference_smiles],
    ["candidate_snapshot_smiles", item.candidate_snapshot_smiles],
    ["candidate_snapshot_runtime", item.candidate_snapshot_runtime],
    ["live_candidate_smiles", item.live_candidate_smiles],
    ["live_candidate_smiles_exact_match", item.live_candidate_smiles_exact_match],
    ["live_matches_candidate_snapshot", item.live_matches_candidate_snapshot],
    ["live_candidate_reason", item.live_candidate_equivalence_reason],
    ["candidate_organic", item.candidate_organic_smiles],
    ["reference_organic", item.reference_organic_smiles],
    ["reference_formula_status", item.reference_formula_check_status],
    ["xyz_formula", item.xyz_formula],
    ["reference_formula_with_h", item.reference_formula_with_h],
    ["reference_formula_mismatch", item.reference_formula_mismatch_detail],
    ["reference_answer_status", item.reference_answer_status],
    ["reference_answer_reason", item.reference_answer_reason],
    ["error", item.error],
  ];
  $("diagnostics").innerHTML = pairs
    .map(([key, value]) => {
      const displayValue = value === null || value === undefined ? "" : String(value);
      return `<dt>${escapeHtml(key)}</dt><dd>${escapeHtml(displayValue)}</dd>`;
    })
    .join("");
}

function statusBadge(status) {
  if (!status) return badge("missing", "bad");
  if (status === "ok") return badge(status, "ok");
  if (status === "not_run") return badge(status);
  return badge(status, "bad");
}

function versionResultCard(label, status, smiles) {
  return `
    <article class="version-card">
      <header>
        <strong>${escapeHtml(label)}</strong>
        ${statusBadge(status)}
      </header>
      <code>${escapeHtml(smiles || "")}</code>
    </article>`;
}

function renderVersionComparison() {
  const item = state.current;
  if (!item) {
    $("versionComparison").innerHTML = "";
    return;
  }
  const py38Status = item.py38_candidate_cpp_status || item.py38_molgr_cpp_status || "";
  const py310Status = item.py310_candidate_cpp_status || item.py310_molgr_cpp_status || "";
  const py38Smiles = item.py38_candidate_cpp_smiles || item.py38_molgr_cpp_smiles || "";
  const py310Smiles = item.py310_candidate_cpp_smiles || item.py310_molgr_cpp_smiles || "";
  const hasCppRows = py38Status || py310Status || py38Smiles || py310Smiles;
  if (!hasCppRows) {
    $("versionComparison").innerHTML = "";
    return;
  }
  const sameStatus = py38Status === py310Status;
  const sameSmiles = py38Smiles && py310Smiles && py38Smiles === py310Smiles;
  const mismatch = item.category === "python_version_mismatch" || !sameStatus || !sameSmiles;
  $("versionComparison").innerHTML = `
    <header>
      <h4>Python 版本重建结果</h4>
      ${badge(mismatch ? "分歧" : "一致", mismatch ? "warn" : "ok")}
    </header>
    <div class="version-grid">
      ${versionResultCard("Python 3.8 · candidate_cpp", py38Status, py38Smiles)}
      ${versionResultCard("Python 3.10 · candidate_cpp", py310Status, py310Smiles)}
    </div>`;
}

async function loadXyz(item) {
  const container = $("viewer3d");
  $("xyzText").textContent = "loading...";
  try {
    const xyz = await fetch(`/api/cases/${encodeURIComponent(item.case_id)}/xyz`).then((r) => {
      if (!r.ok) throw new Error(`XYZ load failed: ${r.status}`);
      return r.text();
    });
    $("xyzText").textContent = xyz;
    render3d(xyz);
  } catch (error) {
    $("xyzText").textContent = error.message;
    container.innerHTML = `<div class="empty">${escapeHtml(error.message)}</div>`;
  }
}

async function loadCandidateSdf(item) {
  const text = $("candidateSdfText");
  const status = $("candidateSdfStatus");
  const modelStatus = $("candidateModelStatus");
  text.textContent = "loading...";
  status.textContent = "";
  modelStatus.textContent = "";
  state.currentCandidateSdf = "";
  state.currentLiveCandidate = null;
  try {
    const data = await api(`/api/cases/${encodeURIComponent(item.case_id)}/candidate-sdf`);
    state.currentLiveCandidate = data;
    item.live_candidate_status = data.live_candidate_status || "";
    item.live_candidate_smiles = data.live_candidate_smiles || "";
    item.live_candidate_smiles_exact_match = data.live_candidate_smiles_exact_match;
    item.live_matches_candidate_snapshot = data.live_matches_candidate_snapshot;
    item.live_candidate_equivalence_method = data.live_candidate_equivalence_method || "";
    item.live_candidate_equivalence_reason = data.live_candidate_equivalence_reason || "";
    renderCaseHeader();
    renderDiagnostics();
    if (!data.available) {
      text.textContent = data.error || "当前代码重建失败";
      status.textContent = "不可用";
      modelStatus.textContent = "不可用";
      renderCandidate3d("");
      return;
    }
    state.currentCandidateSdf = data.sdf || "";
    text.textContent = state.currentCandidateSdf;
    if (data.live_matches_candidate_snapshot === false) {
      status.textContent = "已生成 · 与候选快照不一致";
    } else if (data.live_matches_candidate_snapshot === true) {
      status.textContent = "已生成 · 与候选快照一致";
    } else {
      status.textContent = state.currentCandidateSdf ? "已生成 · 候选快照不可比" : "空结果";
    }
    renderCandidate3d(state.currentCandidateSdf);
    modelStatus.textContent = state.currentCandidateSdf ? "已渲染" : "空结果";
  } catch (error) {
    text.textContent = error.message;
    status.textContent = "错误";
    modelStatus.textContent = "错误";
    renderCandidate3d("");
  }
}

function applyViewerStyle(viewer) {
  if (state.viewStyle === "sphere") {
    viewer.setStyle({}, { sphere: { scale: 0.28 } });
  } else if (state.viewStyle === "line") {
    viewer.setStyle({}, { line: { linewidth: 2 } });
  } else {
    viewer.setStyle({}, { stick: { radius: 0.16 }, sphere: { scale: 0.22 } });
  }
}

function copyViewerPose(sourceViewer, targetViewer) {
  if (
    sourceViewer &&
    targetViewer &&
    typeof sourceViewer.getView === "function" &&
    typeof targetViewer.setView === "function"
  ) {
    targetViewer.setView(sourceViewer.getView());
  }
}

function render3d(xyz) {
  const container = $("viewer3d");
  container.innerHTML = "";
  if (!window.$3Dmol) {
    container.innerHTML = '<div class="empty">3Dmol.js 未加载；仍可查看 XYZ 文本。</div>';
    return;
  }
  requestAnimationFrame(() => {
    const viewer = window.$3Dmol.createViewer(container, { backgroundColor: "white" });
    viewer.addModel(xyz, "xyz");
    applyViewerStyle(viewer);
    viewer.zoomTo();
    if (typeof viewer.resize === "function") {
      viewer.resize();
    }
    viewer.render();
    state.xyzViewer = viewer;
    if (state.currentCandidateSdf) {
      renderCandidate3d(state.currentCandidateSdf);
    }
  });
}

function renderCandidate3d(sdf) {
  const container = $("viewerCandidate3d");
  container.innerHTML = "";
  if (!sdf) {
    container.innerHTML = '<div class="empty">当前 case 没有可展示的候选 3D 拓扑。</div>';
    state.candidateViewer = null;
    return;
  }
  if (!window.$3Dmol) {
    container.innerHTML = '<div class="empty">3Dmol.js 未加载；仍可查看 SDF 文本。</div>';
    state.candidateViewer = null;
    return;
  }
  requestAnimationFrame(() => {
    const viewer = window.$3Dmol.createViewer(container, { backgroundColor: "white" });
    viewer.addModel(sdf, "sdf");
    applyViewerStyle(viewer);
    viewer.zoomTo();
    copyViewerPose(state.xyzViewer, viewer);
    if (typeof viewer.resize === "function") {
      viewer.resize();
    }
    viewer.render();
    state.candidateViewer = viewer;
  });
}

async function loadPair() {
  await Promise.all([
    loadRender(state.primaryKind, "primary"),
    loadRender(state.secondaryKind, "secondary"),
  ]);
}

async function loadRender(kind, slot) {
  const item = state.current;
  const title = slot === "primary" ? $("primaryRenderTitle") : $("secondaryRenderTitle");
  const box = slot === "primary" ? $("primarySvg") : $("secondarySvg");
  const smiles = slot === "primary" ? $("primarySmiles") : $("secondarySmiles");
  title.textContent = labels[kind] || kind;
  box.innerHTML = '<div class="empty">rendering...</div>';
  setImageZoomState(box, title.textContent);
  smiles.textContent = "";
  try {
    const data =
      kind === "candidate"
        ? state.currentLiveCandidate
        : await api(
            `/api/cases/${encodeURIComponent(item.case_id)}/render?kind=${encodeURIComponent(kind)}`,
          );
    const renderError = kind === "candidate" ? data?.render_error || data?.error : data?.error;
    if (!data || renderError) {
      box.innerHTML = `<div class="empty">${escapeHtml(renderError || "当前重建不可用")}</div>`;
    } else {
      box.innerHTML = data.svg || '<div class="empty">empty render</div>';
    }
    setImageZoomState(box, title.textContent);
    smiles.textContent = kind === "candidate" ? data?.live_candidate_smiles || "" : data?.smiles || "";
  } catch (error) {
    box.innerHTML = `<div class="empty">${escapeHtml(error.message)}</div>`;
    setImageZoomState(box, title.textContent);
  }
}

async function saveReview(status) {
  if (!state.current) return;
  $("saveState").textContent = "saving...";
  try {
    if (status === "manual_reference") {
      await readKetcherToForm({ rethrow: true, requireSmiles: true });
    }
    const payload = {
      status,
      corrected_smiles: $("correctedSmiles").value,
      corrected_molblock: $("correctedMolblock").value,
      notes: $("notes").value,
      reviewer: $("reviewer").value,
    };
    localStorage.setItem("moleculeReviewReviewer", payload.reviewer);
    const result = await api(`/api/cases/${encodeURIComponent(state.current.case_id)}/review`, {
      method: "POST",
      body: JSON.stringify(payload),
    });
    $("saveState").textContent = result.fixture
      ? `saved · ${result.fixture.structure_file}`
      : "saved · no fixture";
    await loadStats();
    await loadCase(state.current.case_id);
  } catch (error) {
    $("saveState").textContent = error.message;
  }
}

async function removeCurrentFixture() {
  const item = state.current;
  if (!item?.fixture) return;
  const confirmed = window.confirm(
    `确定移除 ${item.fixture.structure_file}，并将 ${item.case_id} 标记为待复核吗？`,
  );
  if (!confirmed) return;
  await saveReview("needs_followup");
}

function ketcherWindow() {
  const frame = $("ketcherFrame");
  return frame ? frame.contentWindow : null;
}

function getKetcher() {
  const win = ketcherWindow();
  return win && win.ketcher ? win.ketcher : null;
}

async function waitForKetcher(timeoutMs = 20000) {
  const startedAt = Date.now();
  while (Date.now() - startedAt < timeoutMs) {
    const ketcher = getKetcher();
    if (ketcher) {
      state.ketcherLoaded = true;
      $("ketcherStatus").textContent = "Ketcher ready";
      return ketcher;
    }
    await new Promise((resolve) => setTimeout(resolve, 250));
  }
  throw new Error("Ketcher 尚未加载完成");
}

async function setKetcherMolecule(smiles) {
  if (!smiles) {
    $("ketcherStatus").textContent = "当前 case 没有可载入的 SMILES";
    return;
  }
  $("ketcherStatus").textContent = "loading molecule...";
  try {
    const ketcher = await waitForKetcher();
    if (typeof ketcher.setMolecule !== "function") {
      throw new Error("Ketcher API 缺少 setMolecule()");
    }
    await ketcher.setMolecule(smiles);
    $("ketcherStatus").textContent = "molecule loaded";
  } catch (error) {
    $("ketcherStatus").textContent = error.message;
  }
}

async function readKetcherToForm({ rethrow = false, requireSmiles = false } = {}) {
  $("ketcherStatus").textContent = "reading canvas...";
  try {
    const ketcher = await waitForKetcher();
    let molblock = "";
    let molblockError = null;
    if (typeof ketcher.getMolfile === "function") {
      try {
        molblock = String((await ketcher.getMolfile()) || "");
        if (molblock.trim()) $("correctedMolblock").value = molblock;
      } catch (error) {
        molblockError = error;
      }
    }

    let smiles = "";
    if (typeof ketcher.getSmiles === "function") {
      smiles = String((await ketcher.getSmiles()) || "").trim();
      if (smiles) $("correctedSmiles").value = smiles;
    }
    if (requireSmiles && !smiles) throw new Error("人工修图需要 Ketcher SMILES");
    if (!molblock.trim() && !smiles) {
      if (molblockError) throw molblockError;
      throw new Error("Ketcher 画布为空");
    }
    $("ketcherStatus").textContent = "canvas copied to form";
    return { molblock, smiles };
  } catch (error) {
    $("ketcherStatus").textContent = error.message;
    if (rethrow) throw error;
    return { molblock: "", smiles: "" };
  }
}

function currentCandidateOrganicSmiles() {
  if (!state.current) return "";
  return state.current.candidate_organic_smiles || state.current.candidate_smiles || "";
}

function currentReferenceSmiles() {
  if (!state.current) return "";
  return (
    state.current.reference_organic_smiles ||
    state.current.reference_smiles ||
    ""
  );
}

function bindEvents() {
  initializeLayout();
  bindLayoutResizers();
  observeViewerContainer();
  window.addEventListener("load", () => {
    if (window.$3Dmol && $("xyzText").textContent) {
      render3d($("xyzText").textContent);
    }
  });
  $("ketcherFrame").addEventListener("load", () => {
    waitForKetcher(30000)
      .then(() => {
        const workspace = document.querySelector(".workspace");
        if (workspace) workspace.scrollTop = 0;
      })
      .catch((error) => {
        $("ketcherStatus").textContent = error.message;
      });
  });
  $("refreshStats").addEventListener("click", () => {
    loadStats();
    loadCases(true);
  });
  $("categoryFilter").addEventListener("change", () => loadCases(true));
  $("statusFilter").addEventListener("change", () => loadCases(true));
  $("searchBox").addEventListener("input", debounce(() => loadCases(true), 250));
  $("prevPage").addEventListener("click", () => {
    state.offset = Math.max(0, state.offset - state.limit);
    loadCases();
  });
  $("nextPage").addEventListener("click", () => {
    if (state.offset + state.limit < state.total) {
      state.offset += state.limit;
      loadCases();
    }
  });
  $("loadCurrent").addEventListener("click", () => {
    if (state.current) loadCase(state.current.case_id);
  });
  $("openTrace").addEventListener("click", () => {
    if (!state.current) return;
    window.open(
      `/trace/${encodeURIComponent(state.current.case_id)}`,
      "_blank",
      "noopener,noreferrer",
    );
  });
  $("removeFixture").addEventListener("click", removeCurrentFixture);
  $("closeImageLightbox").addEventListener("click", closeImageLightbox);
  $("imageLightbox").addEventListener("click", (event) => {
    if (event.target === $("imageLightbox")) closeImageLightbox();
  });
  document.querySelectorAll(".svg-box").forEach((box) => {
    box.addEventListener("click", () => openImageLightbox(box));
    box.addEventListener("keydown", (event) => {
      if (event.key !== "Enter" && event.key !== " ") return;
      event.preventDefault();
      openImageLightbox(box);
    });
  });
  $("focusKetcher").addEventListener("click", () => {
    $("ketcherPanel").scrollIntoView({ block: "center", behavior: "smooth" });
    $("ketcherFrame").focus();
  });
  $("loadMolgrToKetcher").addEventListener("click", () =>
    setKetcherMolecule(state.currentCandidateSdf || currentCandidateOrganicSmiles()),
  );
  $("loadReferenceToKetcher").addEventListener("click", () => setKetcherMolecule(currentReferenceSmiles()));
  $("readKetcher").addEventListener("click", () => readKetcherToForm());
  document.querySelectorAll(".view-mode").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll(".view-mode").forEach((b) => b.classList.remove("active"));
      button.classList.add("active");
      state.viewStyle = button.dataset.style;
      if ($("xyzText").textContent) render3d($("xyzText").textContent);
      if (state.currentCandidateSdf) renderCandidate3d(state.currentCandidateSdf);
    });
  });
  document.querySelectorAll(".render-kind").forEach((button) => {
    button.addEventListener("click", async () => {
      document.querySelectorAll(".render-kind").forEach((b) => b.classList.remove("active"));
      button.classList.add("active");
      state.primaryKind = button.dataset.kind;
      state.secondaryKind = state.primaryKind === "candidate" ? "reference" : "candidate";
      await loadPair();
    });
  });
  document.querySelectorAll(".decision").forEach((button) => {
    button.addEventListener("click", () => saveReview(button.dataset.status));
  });
}

function debounce(fn, wait) {
  let timer = null;
  return (...args) => {
    clearTimeout(timer);
    timer = setTimeout(() => fn(...args), wait);
  };
}

async function init() {
  bindEvents();
  const requestedCaseId = new URLSearchParams(window.location.search).get("case")?.trim() || "";
  if (requestedCaseId) $("searchBox").value = requestedCaseId;
  await loadStats();
  await loadCases(true);
  if (requestedCaseId) {
    await loadCase(requestedCaseId);
  } else if (state.cases.length) {
    await loadCase(state.cases[0].case_id);
  }
  const workspace = document.querySelector(".workspace");
  if (workspace) workspace.scrollTop = 0;
}

init().catch((error) => {
  document.body.innerHTML = `<pre>${escapeHtml(error.stack || error.message)}</pre>`;
});
