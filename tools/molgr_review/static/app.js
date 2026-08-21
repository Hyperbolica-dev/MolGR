const state = {
  language: localStorage.getItem("moleculeReviewLanguage") === "en" ? "en" : "zh",
  cases: [],
  total: 0,
  offset: 0,
  limit: 80,
  current: null,
  viewStyle: "stick",
  twoDMode: "skeleton",
  caseRequestToken: 0,
  xyzViewer: null,
  candidateViewer: null,
  referenceXyzViewer: null,
  referenceXyzSdf: "",
  referenceXyzFailure: "",
  mappedComparison: {},
  xyzComparisonMode: "raw",
  currentCandidateSdf: "",
  currentLiveCandidate: null,
  ketcherLoaded: false,
  ketcherStatus: null,
  layout: null,
  referenceRenderStatus: "render_failed",
  referenceRenderError: "",
  graphEvidence: {},
  graphEvidenceLoadingCaseId: "",
  xyzLoadStatus: "unknown",
  xyzLoadError: "",
  savingReview: false,
  undoingReview: false,
  reviewHistory: [],
  familyQa: { enabled: false, families: [], progress: {} },
  activeFamilyId: "",
};

const reviewSessionId = (() => {
  const key = "molgrReviewSessionId";
  const storage = globalThis.sessionStorage;
  let value = storage?.getItem(key);
  if (!value) {
    value = globalThis.crypto?.randomUUID?.() || `${Date.now()}-${Math.random()}`;
    storage?.setItem(key, value);
  }
  return value;
})();

const translations = {
  zh: {
    pageTitle: "分子图审核",
    refresh: "刷新",
    languageToggle: "English",
    category: "类别",
    reviewStatus: "审核状态",
    reviewReasonFilter: "审核理由",
    triageBucket: "Triage bucket",
    familyQueue: "Family QA queue",
    calibrationRelation: "校准关系",
    metalStateTransition: "金属状态",
    repeatCount: "重复次数",
    mappingSource: "映射来源",
    transformation: "变换摘要",
    matchesFamily: "符合 family",
    outlierBlocker: "异常 / blocker",
    approveResonance: "批准 · resonance representation",
    approveRedox: "批准 · redox representation",
    rejectSplit: "拒绝批量 / 拆分 family",
    familyQaSafety: "只更新 pending manifest · 不修改正式审核结论。",
    search: "搜索",
    searchPlaceholder: "id 或 row_index",
    all: "全部",
    previousPage: "上一页",
    nextPage: "下一页",
    chooseCase: "选择一个 case",
    currentCase: "当前审核 case",
    openTrace: "打开 Trace",
    reloadCurrent: "重载当前",
    manualConclusion: "人工结论",
    shortcutHint: "快捷键 1–7 · 上一条/下一条 ← →",
    undoReview: "撤回上一条审核",
    recentReviews: "最近审核",
    justReviewed: "刚刚审核: {caseId} → {status}",
    noRecentReviews: "本次 session 尚无审核操作",
    undoComplete: "已撤回 {caseId}",
    triageEvidence: "审核证据",
    mappingConfidence: "mapping",
    fullTraceEvidence: "完整 Trace evidence",
    removeFixture: "移除 fixture",
    removeFixtureTitle: "移除 {file} 并将 {caseId} 标记为待复核",
    reviewer: "审核理由",
    candidateGraph: "候选图为准",
    referenceGraph: "参考图为准",
    acceptBoth: "接受两者",
    manualReference: "人工修图为准",
    referenceAnswerWrong: "参考答案错误",
    needsFollowup: "待复核",
    skip: "跳过",
    xyz3d: "3D XYZ",
    referenceXyz: "参考 XYZ",
    referenceXyzHint: "Reference graph on source XYZ",
    referenceXyzUnavailable: "Reference XYZ unavailable",
    referenceXyzUnreliable: "原子对应关系不可靠。",
    referenceXyzAmbiguous: "原子对应关系存在歧义。",
    representativeMappingWarning: "代表性映射 — 原子对应关系不唯一。",
    ambiguityType: "歧义类型",
    unavailableReason: "不可用原因",
    ambiguityLocation: "歧义位置",
    ambiguityAlternatives: "可能对应",
    ambiguityLocationUnknown: "无法可靠定位；映射枚举在确认唯一对应关系前已停止。",
    mappingAlternative: "映射 {index}",
    mappingAmbiguityTypes: {
      multiple_valid_mappings: "多个有效映射",
      symmetry_equivalent_atoms: "对称等价原子",
      mapping_enumeration_truncated: "映射枚举被截断",
      mapping_timeout: "映射计算超时",
      ambiguous_mapping: "原子映射不唯一",
    },
    mappingAmbiguityReasons: {
      multiple_equally_valid_atom_mappings: "存在多个同等有效的原子映射，因此无法构造唯一的 Reference XYZ。",
      symmetry_equivalent_atoms_prevent_unique_correspondence: "对称等价原子导致对应关系不唯一，因此无法构造唯一的 Reference XYZ。",
      mapping_enumeration_truncated_before_unique_correspondence: "在确认唯一对应关系之前映射枚举已被截断。",
      mapping_timeout_before_unique_correspondence: "映射计算在确认唯一对应关系之前超时。",
      unique_atom_correspondence_not_established: "无法建立唯一的 Candidate–Reference 原子对应关系。",
    },
    referenceXyzMissing: "Reference graph 缺失。",
    rawXyz: "Raw XYZ",
    mappedComparison: "Mapped comparison",
    mappedDonorPreserved: "mapped donor preserved",
    donorSwappedMappedLigand: "donor swapped within mapped ligand",
    skeletonMode: "骨架",
    hydrogenMode: "H 分配",
    inputXyz: "输入 XYZ",
    currentCandidateTopology3d: "当前候选重建 3D 拓扑",
    currentCandidateSdf: "当前候选重建 SDF",
    xyzText: "XYZ 文本",
    topologyCompare: "拓扑对比",
    currentCandidate: "当前候选",
    reference: "参考图",
    candidateOrganic: "候选 Organic graph",
    referenceOrganic: "参考 Organic graph",
    manualEdit: "人工修图",
    ketcherCorrection: "Ketcher 修正拓扑",
    loadCandidateSdf: "载入候选 SDF",
    loadReference: "载入参考图",
    readCanvas: "读取画布",
    correctedSmiles: "修正 SMILES",
    correctedMolblock: "修正 Molfile / Molblock",
    notes: "备注 / 证据",
    reviewerSummary: "Reviewer Summary / 审核关注项",
    datasetRuntime: "数据集 / 运行环境",
    reviewStructures: "审核结构",
    reviewDetailsJump: "审核详情",
    reviewStructuresHint: "优先核对输入几何、当前候选图与参考图。",
    xyzGeometry: "原始 3D 几何",
    additionalViews: "其他结构视图",
    reviewerDetails: "Reviewer details / 审核详情",
    smilesGraphSummary: "SMILES & graph summary",
    atomBondInspector: "Atom / bond inspector",
    provenanceSnapshot: "Provenance / snapshot",
    showFullGraph: "查看完整 graph",
    metric: "指标",
    element: "元素",
    formalCharge: "形式电荷",
    radicalElectrons: "自由基电子",
    explicitH: "显式 H",
    implicitH: "隐式 H",
    neighbours: "邻接原子",
    bond: "键",
    bondType: "类型",
    currentCandidateSmiles: "Current Candidate SMILES",
    candidateSnapshotSmiles: "Candidate snapshot SMILES",
    referenceSmiles: "Reference SMILES",
    totalFormalCharge: "总形式电荷",
    atomCount: "原子数",
    explicitHCount: "显式 H 数",
    totalRadicalElectrons: "总自由基电子数",
    metals: "金属及形式电荷",
    graphLoading: "正在读取 graph…",
    graphUnavailable: "Graph unavailable",
    noRelevantAtoms: "没有符合默认筛选的原子",
    currentCandidateStatus: "当前候选状态",
    currentVsSnapshot: "Current vs snapshot",
    equivalenceMethod: "等价性方法",
    equivalenceReason: "等价性说明",
    snapshotRuntime: "Snapshot runtime",
    qTooltip: "分子总电荷",
    multiplicityTooltip: "自旋多重度",
    radicalsTooltip: "显式自由基电子总数",
    formulaTooltip: "XYZ 与 graph 分子式一致性",
    electronStateRoundtripNote: "LONE_PAIR_COUNT_PROP 丢失可改变价态语义",
    caseDetails: "Developer details / 开发信息",
    assessmentDetails: "Assessment / Reference 详情",
    runtimeDetails: "Benchmark / Runtime",
    otherMetadata: "Developer / raw metadata",
    notProvided: "未提供",
    present: "存在",
    missing: "缺失",
    invalid: "无效",
    valid: "有效",
    unknownValidity: "存在；有效性未提供",
    currentDecision: "当前已保存：{status}",
    noSavedDecision: "当前未保存审核结论",
    reviewedBy: "审核理由：{reviewer}",
    updatedAt: "更新于：{updatedAt}",
    referenceMissingNotice: "当前 payload 未提供 reference SMILES。",
    referenceInvalidNotice: "当前 reference SMILES 无法解析。",
    referenceRenderFailed: "Reference render failed",
    focusCandidateFailure: "重点检查：当前候选重建不可用。",
    focusReference: "重点检查：Reference 缺失、无效或渲染失败，优先核对 XYZ 与 Candidate。",
    focusReferenceComparison: "重点检查：Reference 比较未完成；这不是分子式错误。",
    focusFormula: "重点检查：XYZ 与 Reference 的分子式状态异常。",
    focusAssessment: "重点检查：当前 case 的 assessability 受限。",
    focusSnapshot: "重点检查：当前重建与 candidate snapshot 不一致。",
    focusComparison: "重点检查：对照 XYZ、Candidate 与 Reference 的成键和电荷。",
    summaryReference: "Reference",
    summaryCandidate: "Candidate",
    summaryCharge: "Charge",
    summaryMultiplicity: "Multiplicity",
    summaryRadicals: "Radicals",
    summaryFormula: "Formula",
    summaryAssessability: "Assessability",
    summarySnapshot: "Snapshot vs current",
    summaryFixture: "Fixture",
    xyzUnavailable: "XYZ 无法加载",
    candidateUnavailableBecauseXyz: "当前候选因此不可用",
    candidateUnavailable: "当前候选不可用",
    technicalDetail: "技术信息：{detail}",
    statusOk: "✓",
    statusMissingCompact: "缺失",
    statusInvalidCompact: "无效",
    statusUnavailableCompact: "不可用",
    statusUnknownCompact: "—",
    formulaOk: "✓",
    snapshotCurrent: "current",
    snapshotDifferent: "不同",
    snapshotUnknown: "—",
    queueLabel: "Queue",
    currentCandidateLabel: "Current",
    currentVsSnapshotLabel: "Current vs snapshot",
    availableCompact: "available",
    missingCompact: "missing",
    invalidCompact: "invalid",
    renderFailedCompact: "render failed",
    unavailableCompact: "unavailable",
    sameCompact: "same",
    differentCompact: "different",
    formulaLabel: "formula",
    diagnosticLabel: "诊断",
    diagnosticReasonLabel: "原因",
    referenceDiagnostics: {
      missing_reference: "缺少 Reference",
      equivalence_timeout: "等价性比较超时",
      candidate_reparse_failure: "Candidate 重新解析失败",
      formula_mismatch: "分子式不匹配",
      comparison_skipped: "比较未完成",
    },
    structure: "结构图",
    close: "关闭",
    closeZoomedImage: "关闭放大图片",
    resizeSidebar: "调整左侧队列宽度",
    resizeCompare: "调整 3D 与对比栏宽度",
    resizeKetcher: "调整 Ketcher 高度",
    cases: "cases",
    row: "行",
    unreviewed: "未审核",
    noFixture: "无 fixture",
    fixture: "fixture · {kind}",
    candidateSnapshot: "候选快照 · {status}",
    candidateSnapshotMismatch: "当前重建 != 候选快照",
    traceTitle: "在新窗口打开 {caseId} 的 Trace",
    selectCaseFirst: "请先选择 case",
    statsAll: "全部",
    statsUnreviewed: "未审核",
    statsCandidate: "候选图为准",
    statsReference: "参考图为准",
    statsBoth: "接受两者",
    statsManual: "人工修图",
    statsWrong: "参考答案错误",
    statusMissing: "missing",
    statusLoading: "加载中...",
    unavailable: "不可用",
    generatedSnapshotMismatch: "已生成 · 与候选快照不一致",
    generatedSnapshotMatch: "已生成 · 与候选快照一致",
    generatedSnapshotIncomparable: "已生成 · 候选快照不可比",
    emptyResult: "空结果",
    rendered: "已渲染",
    error: "错误",
    rendering: "渲染中...",
    emptyRender: "空渲染",
    reconstructionUnavailable: "当前重建不可用",
    emptyCandidate3d: "当前 case 没有可展示的候选 3D 拓扑。",
    threeDmolUnavailableXyz: "3Dmol.js 未加载；仍可查看 XYZ 文本。",
    threeDmolUnavailableSdf: "3Dmol.js 未加载；仍可查看 SDF 文本。",
    saving: "保存中...",
    savedFixture: "已保存 · {file}",
    savedNoFixture: "已保存 · 无 fixture",
    ketcherReady: "Ketcher 已就绪",
    noLoadableSmiles: "当前 case 没有可载入的 SMILES",
    loadingMolecule: "加载分子中...",
    moleculeLoaded: "分子已载入",
    readingCanvas: "读取画布中...",
    canvasCopied: "画布内容已复制到表单",
    ketcherNotReady: "Ketcher 尚未加载完成",
    ketcherNoSetMolecule: "Ketcher API 缺少 setMolecule()",
    manualSmilesRequired: "人工修图需要 Ketcher SMILES",
    canvasEmpty: "Ketcher 画布为空",
    removeFixtureConfirm: "确定移除 {file}，并将 {caseId} 标记为待复核吗？",
    zoomImage: "放大{label}",
    versionComparison: "Python 版本重建结果",
    disagreement: "分歧",
    consistent: "一致",
    diagnostics: {
      total_charge: "总电荷",
      spin_multiplicity: "多重度",
      total_radical_electrons: "总自由基电子数",
      reference_smiles: "参考 SMILES",
      candidate_smiles: "候选 SMILES",
      candidate_snapshot_smiles: "候选快照 SMILES",
      candidate_snapshot_runtime: "候选快照运行时",
      live_candidate_smiles: "当前候选 SMILES",
      live_candidate_smiles_exact_match: "当前候选 SMILES 精确匹配",
      live_matches_candidate_snapshot: "当前候选匹配候选快照",
      live_candidate_reason: "当前候选等价性原因",
      candidate_organic: "候选图 organic",
      reference_organic: "参考图 organic",
      reference_formula_status: "参考分子式状态",
      reference_diagnostic_group: "Reference 诊断分类",
      reference_diagnostic_reason: "Reference 诊断原因",
      xyz_formula: "XYZ 分子式",
      reference_formula_with_h: "参考分子式（含氢）",
      reference_formula_mismatch: "参考分子式差异",
      reference_answer_status: "参考答案状态",
      reference_answer_reason: "参考答案原因",
      accuracy_assessment_status: "准确性评估状态",
      accuracy_assessment_reason: "准确性评估原因",
      tmqmg_answer_assessment: "tmQMg 答案评估",
      molgr_answer_assessment: "MolGR 答案评估",
      error: "错误",
    },
    categories: {
      graph_not_equivalent: "图不等价",
      missing_reference_smiles: "缺参考",
      candidate_failed: "候选生成失败",
      backend_mismatch: "后端分歧",
      python_version_mismatch: "版本分歧",
      reference_not_comparable: "参考不可比",
      reference_formula_mismatch: "参考氢数不守恒",
      no_clear_evidence_boron_cluster: "硼簇结构不可判定",
    },
    status: {
      accept_candidate: "候选图为准",
      accept_reference: "参考图为准",
      accept_both: "接受两者",
      manual_reference: "人工修图为准",
      reference_answer_wrong: "参考答案错误",
      needs_followup: "待复核",
      skip: "跳过",
      unreviewed: "未审核",
    },
  },
  en: {
    pageTitle: "Molecule Graph Review",
    refresh: "Refresh",
    languageToggle: "中文",
    category: "Category",
    reviewStatus: "Review status",
    reviewReasonFilter: "Review reason",
    triageBucket: "Triage bucket",
    familyQueue: "Family QA queue",
    calibrationRelation: "Calibration relation",
    metalStateTransition: "Metal state",
    repeatCount: "Repeat count",
    mappingSource: "Mapping source",
    transformation: "Transformation",
    matchesFamily: "Matches family",
    outlierBlocker: "Outlier / blocker",
    approveResonance: "Approve · resonance representation",
    approveRedox: "Approve · redox representation",
    rejectSplit: "Reject batch / split family",
    familyQaSafety: "Pending manifest only · authoritative reviews are unchanged.",
    search: "Search",
    searchPlaceholder: "id or row_index",
    all: "All",
    previousPage: "Previous",
    nextPage: "Next",
    chooseCase: "Select a case",
    currentCase: "Current review case",
    openTrace: "Open Trace",
    reloadCurrent: "Reload current",
    manualConclusion: "Manual decision",
    shortcutHint: "Shortcuts 1–7 · previous/next ← →",
    undoReview: "Undo last review",
    recentReviews: "Recent reviews",
    justReviewed: "Just reviewed: {caseId} → {status}",
    noRecentReviews: "No reviews in this session",
    undoComplete: "Undid {caseId}",
    triageEvidence: "Review evidence",
    mappingConfidence: "mapping",
    fullTraceEvidence: "Full Trace evidence",
    removeFixture: "Remove fixture",
    removeFixtureTitle: "Remove {file} and mark {caseId} for follow-up",
    reviewer: "Review reason",
    candidateGraph: "Accept candidate",
    referenceGraph: "Accept reference",
    acceptBoth: "Accept both",
    manualReference: "Use manual edit",
    referenceAnswerWrong: "Reference answer wrong",
    needsFollowup: "Needs follow-up",
    skip: "Skip",
    xyz3d: "3D XYZ",
    referenceXyz: "Reference XYZ",
    referenceXyzHint: "Reference graph on source XYZ",
    referenceXyzUnavailable: "Reference XYZ unavailable",
    referenceXyzUnreliable: "Atom correspondence is not reliable.",
    referenceXyzAmbiguous: "Atom correspondence is ambiguous.",
    representativeMappingWarning: "Representative mapping — atom correspondence is not unique.",
    ambiguityType: "Ambiguity type",
    unavailableReason: "Reason",
    ambiguityLocation: "Ambiguous location",
    ambiguityAlternatives: "Possible correspondences",
    ambiguityLocationUnknown: "Not reliably localized; mapping search stopped before a unique correspondence was established.",
    mappingAlternative: "Mapping {index}",
    mappingAmbiguityTypes: {
      multiple_valid_mappings: "multiple valid mappings",
      symmetry_equivalent_atoms: "symmetry-equivalent atoms",
      mapping_enumeration_truncated: "mapping enumeration truncated",
      mapping_timeout: "mapping timeout",
      ambiguous_mapping: "ambiguous atom mapping",
    },
    mappingAmbiguityReasons: {
      multiple_equally_valid_atom_mappings: "Multiple equally valid atom mappings exist, so a unique Reference XYZ cannot be constructed.",
      symmetry_equivalent_atoms_prevent_unique_correspondence: "Symmetry-equivalent atoms prevent a unique correspondence, so a unique Reference XYZ cannot be constructed.",
      mapping_enumeration_truncated_before_unique_correspondence: "Mapping enumeration was truncated before a unique correspondence could be established.",
      mapping_timeout_before_unique_correspondence: "Mapping timed out before a unique correspondence could be established.",
      unique_atom_correspondence_not_established: "A unique Candidate–Reference atom correspondence could not be established.",
    },
    referenceXyzMissing: "Reference graph is missing.",
    rawXyz: "Raw XYZ",
    mappedComparison: "Mapped comparison",
    mappedDonorPreserved: "mapped donor preserved",
    donorSwappedMappedLigand: "donor swapped within mapped ligand",
    skeletonMode: "Skeleton",
    hydrogenMode: "Show H",
    inputXyz: "Input XYZ",
    currentCandidateTopology3d: "Current candidate 3D topology",
    currentCandidateSdf: "Current candidate SDF",
    xyzText: "XYZ text",
    topologyCompare: "Topology comparison",
    currentCandidate: "Current candidate",
    reference: "Reference",
    candidateOrganic: "Candidate Organic graph",
    referenceOrganic: "Reference Organic graph",
    manualEdit: "Manual editing",
    ketcherCorrection: "Ketcher topology correction",
    loadCandidateSdf: "Load candidate SDF",
    loadReference: "Load reference",
    readCanvas: "Read canvas",
    correctedSmiles: "Corrected SMILES",
    correctedMolblock: "Corrected Molfile / Molblock",
    notes: "Notes / evidence",
    reviewerSummary: "Reviewer Summary",
    datasetRuntime: "Dataset / runtime",
    reviewStructures: "Review structures",
    reviewDetailsJump: "Review details",
    reviewStructuresHint: "Check the input geometry, current candidate, and reference first.",
    xyzGeometry: "Original 3D geometry",
    additionalViews: "Additional structure views",
    reviewerDetails: "Reviewer details",
    smilesGraphSummary: "SMILES & graph summary",
    atomBondInspector: "Atom / bond inspector",
    provenanceSnapshot: "Provenance / snapshot",
    showFullGraph: "Show full graph",
    metric: "Metric",
    element: "Element",
    formalCharge: "Formal charge",
    radicalElectrons: "Radical e−",
    explicitH: "Explicit H",
    implicitH: "Implicit H",
    neighbours: "Neighbours",
    bond: "Bond",
    bondType: "Type",
    currentCandidateSmiles: "Current Candidate SMILES",
    candidateSnapshotSmiles: "Candidate snapshot SMILES",
    referenceSmiles: "Reference SMILES",
    totalFormalCharge: "Total formal charge",
    atomCount: "Atom count",
    explicitHCount: "Explicit H count",
    totalRadicalElectrons: "Total radical electrons",
    metals: "Metals and formal charges",
    graphLoading: "Loading graph…",
    graphUnavailable: "Graph unavailable",
    noRelevantAtoms: "No atoms match the default filter",
    currentCandidateStatus: "Current candidate status",
    currentVsSnapshot: "Current vs snapshot",
    equivalenceMethod: "Equivalence method",
    equivalenceReason: "Equivalence reason",
    snapshotRuntime: "Snapshot runtime",
    qTooltip: "Total molecular charge",
    multiplicityTooltip: "Spin multiplicity",
    radicalsTooltip: "Total explicit radical electrons",
    formulaTooltip: "XYZ / graph molecular formula consistency",
    electronStateRoundtripNote: "LONE_PAIR_COUNT_PROP loss can change valence semantics",
    caseDetails: "Developer details",
    assessmentDetails: "Assessment / Reference details",
    runtimeDetails: "Benchmark / Runtime",
    otherMetadata: "Developer / raw metadata",
    notProvided: "Not provided",
    present: "Present",
    missing: "Missing",
    invalid: "Invalid",
    valid: "Valid",
    unknownValidity: "Present; validity not provided",
    currentDecision: "Saved decision: {status}",
    noSavedDecision: "No review decision saved",
    reviewedBy: "Review reason: {reviewer}",
    updatedAt: "Updated: {updatedAt}",
    referenceMissingNotice: "The current payload does not provide a reference SMILES.",
    referenceInvalidNotice: "The reference SMILES could not be parsed.",
    referenceRenderFailed: "Reference render failed",
    focusCandidateFailure: "Focus: the current candidate reconstruction is unavailable.",
    focusReference: "Focus: reference is missing, invalid, or failed to render; compare XYZ and Candidate first.",
    focusReferenceComparison: "Focus: Reference comparison did not complete; this is not a formula error.",
    focusFormula: "Focus: the XYZ and Reference formula status is abnormal.",
    focusAssessment: "Focus: assessability is limited for this case.",
    focusSnapshot: "Focus: current reconstruction differs from the candidate snapshot.",
    focusComparison: "Focus: compare bonding and charge across XYZ, Candidate, and Reference.",
    summaryReference: "Reference",
    summaryCandidate: "Candidate",
    summaryCharge: "Charge",
    summaryMultiplicity: "Multiplicity",
    summaryRadicals: "Radicals",
    summaryFormula: "Formula",
    summaryAssessability: "Assessability",
    summarySnapshot: "Snapshot vs current",
    summaryFixture: "Fixture",
    xyzUnavailable: "XYZ could not be loaded",
    candidateUnavailableBecauseXyz: "The current candidate is therefore unavailable",
    candidateUnavailable: "Current candidate unavailable",
    technicalDetail: "Technical detail: {detail}",
    statusOk: "✓",
    statusMissingCompact: "missing",
    statusInvalidCompact: "invalid",
    statusUnavailableCompact: "unavailable",
    statusUnknownCompact: "—",
    formulaOk: "✓",
    snapshotCurrent: "current",
    snapshotDifferent: "different",
    snapshotUnknown: "—",
    queueLabel: "Queue",
    currentCandidateLabel: "Current",
    currentVsSnapshotLabel: "Current vs snapshot",
    availableCompact: "available",
    missingCompact: "missing",
    invalidCompact: "invalid",
    renderFailedCompact: "render failed",
    unavailableCompact: "unavailable",
    sameCompact: "same",
    differentCompact: "different",
    formulaLabel: "formula",
    diagnosticLabel: "Diagnosis",
    diagnosticReasonLabel: "Reason",
    referenceDiagnostics: {
      missing_reference: "Missing reference",
      equivalence_timeout: "Equivalence comparison timeout",
      candidate_reparse_failure: "Candidate reparse failure",
      formula_mismatch: "Formula mismatch",
      comparison_skipped: "Comparison incomplete",
    },
    structure: "Structure",
    close: "Close",
    closeZoomedImage: "Close enlarged image",
    resizeSidebar: "Resize case queue",
    resizeCompare: "Resize 3D and comparison panes",
    resizeKetcher: "Resize Ketcher height",
    cases: "cases",
    row: "row",
    unreviewed: "Unreviewed",
    noFixture: "No fixture",
    fixture: "fixture · {kind}",
    candidateSnapshot: "Candidate snapshot · {status}",
    candidateSnapshotMismatch: "Live reconstruction != candidate snapshot",
    traceTitle: "Open Trace for {caseId} in a new window",
    selectCaseFirst: "Select a case first",
    statsAll: "All",
    statsUnreviewed: "Unreviewed",
    statsCandidate: "Accept candidate",
    statsReference: "Accept reference",
    statsBoth: "Accept both",
    statsManual: "Manual edit",
    statsWrong: "Reference answer wrong",
    statusMissing: "missing",
    statusLoading: "Loading...",
    unavailable: "Unavailable",
    generatedSnapshotMismatch: "Generated · differs from candidate snapshot",
    generatedSnapshotMatch: "Generated · matches candidate snapshot",
    generatedSnapshotIncomparable: "Generated · candidate snapshot incomparable",
    emptyResult: "Empty result",
    rendered: "Rendered",
    error: "Error",
    rendering: "Rendering...",
    emptyRender: "Empty render",
    reconstructionUnavailable: "Current reconstruction unavailable",
    emptyCandidate3d: "No candidate 3D topology is available for this case.",
    threeDmolUnavailableXyz: "3Dmol.js is not loaded; XYZ text remains available.",
    threeDmolUnavailableSdf: "3Dmol.js is not loaded; SDF text remains available.",
    saving: "Saving...",
    savedFixture: "Saved · {file}",
    savedNoFixture: "Saved · no fixture",
    ketcherReady: "Ketcher ready",
    noLoadableSmiles: "This case has no loadable SMILES",
    loadingMolecule: "Loading molecule...",
    moleculeLoaded: "Molecule loaded",
    readingCanvas: "Reading canvas...",
    canvasCopied: "Canvas copied to form",
    ketcherNotReady: "Ketcher has not finished loading",
    ketcherNoSetMolecule: "Ketcher API is missing setMolecule()",
    manualSmilesRequired: "Manual editing requires a Ketcher SMILES",
    canvasEmpty: "Ketcher canvas is empty",
    removeFixtureConfirm: "Remove {file} and mark {caseId} for follow-up?",
    zoomImage: "Enlarge {label}",
    versionComparison: "Python version reconstruction",
    disagreement: "Mismatch",
    consistent: "Consistent",
    diagnostics: {
      total_charge: "Total charge",
      spin_multiplicity: "Spin multiplicity",
      total_radical_electrons: "Total radical electrons",
      reference_smiles: "Reference SMILES",
      candidate_smiles: "Candidate SMILES",
      candidate_snapshot_smiles: "Candidate snapshot SMILES",
      candidate_snapshot_runtime: "Candidate snapshot runtime",
      live_candidate_smiles: "Live candidate SMILES",
      live_candidate_smiles_exact_match: "Live candidate exact SMILES match",
      live_matches_candidate_snapshot: "Live candidate matches snapshot",
      live_candidate_reason: "Live candidate equivalence reason",
      candidate_organic: "Candidate organic",
      reference_organic: "Reference organic",
      reference_formula_status: "Reference formula status",
      reference_diagnostic_group: "Reference diagnostic group",
      reference_diagnostic_reason: "Reference diagnostic reason",
      xyz_formula: "XYZ formula",
      reference_formula_with_h: "Reference formula with H",
      reference_formula_mismatch: "Reference formula mismatch",
      reference_answer_status: "Reference answer status",
      reference_answer_reason: "Reference answer reason",
      accuracy_assessment_status: "Accuracy assessment status",
      accuracy_assessment_reason: "Accuracy assessment reason",
      tmqmg_answer_assessment: "tmQMg answer assessment",
      molgr_answer_assessment: "MolGR answer assessment",
      error: "Error",
    },
    categories: {
      graph_not_equivalent: "Graph not equivalent",
      missing_reference_smiles: "Missing reference",
      candidate_failed: "Candidate reconstruction failed",
      backend_mismatch: "Backend mismatch",
      python_version_mismatch: "Python version mismatch",
      reference_not_comparable: "Reference not comparable",
      reference_formula_mismatch: "Reference H-count mismatch",
      no_clear_evidence_boron_cluster: "Boron cluster not assessable",
    },
    status: {
      accept_candidate: "Accept candidate",
      accept_reference: "Accept reference",
      accept_both: "Accept both",
      manual_reference: "Use manual edit",
      reference_answer_wrong: "Reference answer wrong",
      needs_followup: "Needs follow-up",
      skip: "Skip",
      unreviewed: "Unreviewed",
    },
  },
};

function tr(key, fallback = key) {
  const dictionary = translations[state.language] || translations.zh;
  const value = key.split(".").reduce((current, part) => current?.[part], dictionary);
  return value === undefined ? fallback : value;
}

function msg(key, values = {}) {
  return Object.entries(values).reduce(
    (text, [name, value]) => text.replaceAll(`{${name}}`, String(value)),
    tr(key),
  );
}

function categoryLabel(category) {
  return tr(`categories.${category}`, category);
}

function statusLabel(status) {
  return status ? tr(`status.${status}`, status) : tr("status.unreviewed");
}

function applyLanguage() {
  document.documentElement.lang = state.language === "en" ? "en" : "zh-CN";
  document.title = tr("pageTitle");
  document.querySelectorAll("[data-i18n]").forEach((element) => {
    element.textContent = tr(element.dataset.i18n);
  });
  document.querySelectorAll("[data-i18n-placeholder]").forEach((element) => {
    element.placeholder = tr(element.dataset.i18nPlaceholder);
  });
  document.querySelectorAll("[data-i18n-title]").forEach((element) => {
    element.title = tr(element.dataset.i18nTitle);
  });
  document.querySelectorAll("[data-i18n-aria-label]").forEach((element) => {
    element.setAttribute("aria-label", tr(element.dataset.i18nAriaLabel));
  });
  localizeTriageFilterOptions();
  renderCaseList();
  renderMappedComparisonNote();
  if (state.current) {
    renderCaseHeader();
    renderReviewerSummary();
    renderVersionComparison();
    renderDiagnostics();
    renderReviewerDetails();
    renderCandidateSdfStatus();
    loadPair(state.caseRequestToken);
  }
  renderKetcherStatus();
}

function toggleLanguage() {
  state.language = state.language === "zh" ? "en" : "zh";
  localStorage.setItem("moleculeReviewLanguage", state.language);
  applyLanguage();
  loadStats().catch((error) => {
    console.error("Failed to localize sidebar counts", error);
  });
}

function renderKetcherStatus() {
  const status = state.ketcherStatus;
  if (!status) return;
  $("ketcherStatus").textContent = status.key ? msg(status.key, status.values) : status.text;
}

function setKetcherStatus(key, values = {}, text = "") {
  state.ketcherStatus = { key, values, text };
  renderKetcherStatus();
}

function localizedError(key) {
  const error = new Error(tr(key));
  error.i18nKey = key;
  return error;
}

const labels = {
  candidate: "currentCandidate",
  reference: "reference",
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
    [state.xyzViewer, state.referenceXyzViewer, state.candidateViewer].forEach((viewer) => {
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
  const containers = [$("viewer3d"), $("referenceViewer3d"), $("viewerCandidate3d")].filter(Boolean);
  if (!containers.length) return;
  viewerResizeObserver = new ResizeObserver(() => queueViewerResize());
  containers.forEach((container) => viewerResizeObserver.observe(container));
}

async function api(path, options = {}) {
  const response = await fetch(path, {
    ...options,
    headers: {
      "Content-Type": "application/json",
      "X-MolGR-Review-Session": reviewSessionId,
      ...(options.headers || {}),
    },
  });
  const contentType = response.headers.get("content-type") || "";
  const payload = contentType.includes("application/json") ? await response.json() : await response.text();
  if (!response.ok) {
    const message = typeof payload === "string" ? payload : payload.error || response.statusText;
    const error = new Error(message);
    error.httpStatus = response.status;
    error.payloadError = typeof payload === "object" && payload ? payload.error || "" : "";
    throw error;
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
  box.setAttribute("aria-label", msg("zoomImage", { label }));
  box.dataset.zoomLabel = label;
}

function openImageLightbox(box) {
  const image = box.querySelector("svg, img");
  if (!image) return;
  const dialog = $("imageLightbox");
  $("imageLightboxTitle").textContent = box.dataset.zoomLabel || tr("structure");
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
  if (category === "no_clear_evidence_boron_cluster") return "warn";
  if (category === "candidate_failed") return "bad";
  return "";
}

function reviewStatusKind(status) {
  return status === "reference_answer_wrong" ? "bad" : "ok";
}

async function loadStats() {
  const stats = await api("/api/stats");
  const categoryEntries = Object.entries(stats.categories || {})
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([key, value]) => [categoryLabel(key), value]);
  const entries = [
    [tr("statsAll"), Object.values(stats.categories || {}).reduce((a, b) => a + b, 0)],
    [tr("statsUnreviewed"), (stats.review_statuses || {}).unreviewed || 0],
    [tr("statsCandidate"), (stats.review_statuses || {}).accept_candidate || 0],
    [tr("statsReference"), (stats.review_statuses || {}).accept_reference || 0],
    [tr("statsBoth"), (stats.review_statuses || {}).accept_both || 0],
    [tr("statsManual"), (stats.review_statuses || {}).manual_reference || 0],
    [tr("statsWrong"), (stats.review_statuses || {}).reference_answer_wrong || 0],
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
  renderTriageBucketOptions(stats.triage_buckets || {});
}

function renderTriageBucketOptions(triageBuckets) {
  const triageFilter = $("triageFilter");
  const selectedBucket = triageFilter.value;
  const triageEntries = Object.entries(triageBuckets).sort(([left], [right]) =>
    left.localeCompare(right),
  );
  if (selectedBucket && !triageEntries.some(([bucket]) => bucket === selectedBucket)) {
    triageEntries.push([selectedBucket, 0]);
  }
  $("triageFilterField").hidden = triageEntries.length === 0;
  triageFilter.innerHTML = `<option value="">${escapeHtml(tr("all"))}</option>${triageEntries
    .map(
      ([bucket, count]) =>
        `<option value="${escapeHtml(bucket)}" data-count="${count}">${escapeHtml(triageBucketLabel(bucket))} (${count})</option>`,
    )
    .join("")}`;
  if (triageEntries.some(([bucket]) => bucket === selectedBucket)) {
    triageFilter.value = selectedBucket;
  }
}

async function loadReviewReasons() {
  const data = await api("/api/review-reasons");
  const items = Array.isArray(data.items) ? data.items : [];
  const reasonFilter = $("reviewReasonFilter");
  const selectedReason = reasonFilter.value;
  $("reviewReasonOptions").innerHTML = items
    .filter((item) => hasValue(item.reviewer))
    .map((item) => {
      const reviewer = String(item.reviewer).trim();
      const count = Number(item.count) || 0;
      return `<option value="${escapeHtml(reviewer)}">${escapeHtml(`${reviewer} (${count})`)}</option>`;
    })
    .join("");
  reasonFilter.innerHTML = `<option value="">${escapeHtml(tr("all"))}</option>${items
    .filter((item) => hasValue(item.reviewer))
    .map((item) => {
      const reviewer = String(item.reviewer).trim();
      const count = Number(item.count) || 0;
      return `<option value="${escapeHtml(reviewer)}">${escapeHtml(`${reviewer} (${count})`)}</option>`;
    })
    .join("")}`;
  if (items.some((item) => String(item.reviewer).trim() === selectedReason)) {
    reasonFilter.value = selectedReason;
  }
}

function activeFamily() {
  return state.familyQa.families.find((family) => family.family_id === state.activeFamilyId) || null;
}

function familyOptionLabel(family) {
  return `${family.family_id} (${family.family_size} cases / ${family.representatives.length} reps)`;
}

function signedNumber(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value ?? "—");
  return number > 0 ? `+${number}` : String(number).replace("-", "−");
}

function bondSymbol(type) {
  return { single: "–", double: "=", triple: "≡", aromatic: ":", dative: "→" }[type] || `(${type})`;
}

function transformationSummary(raw, repeatCount = 1) {
  let value;
  try {
    value = typeof raw === "string" ? JSON.parse(raw) : raw;
  } catch (error) {
    return tr("unknown");
  }
  if (!value || typeof value !== "object") return tr("unknown");
  const lines = [];
  const repeats = Math.max(1, Number(repeatCount) || 1);
  if (hasValue(value.metal_charge_delta) && Number(value.metal_charge_delta) !== 0) {
    lines.push(`metal Δq ${signedNumber(value.metal_charge_delta)}`);
  }
  (value.charge_transitions_per_unit || []).forEach(([element, before, after, count]) => {
    const total = (Number(count) || 1) * repeats;
    lines.push(`${total > 1 ? `${total} × ` : ""}${element} ${signedNumber(before)} → ${signedNumber(after)}`);
  });
  (value.bond_transitions_per_unit || []).forEach(([atoms, before, after, count]) => {
    const [left, right] = atoms || ["?", "?"];
    const total = (Number(count) || 1) * repeats;
    lines.push(`${total > 1 ? `${total} × ` : ""}${left}${bondSymbol(before)}${right} → ${left}${bondSymbol(after)}${right}`);
  });
  if (hasValue(value.ligand_charge_compensation)) {
    lines.push(`ligand charge compensation = ${signedNumber(value.ligand_charge_compensation)}`);
  }
  return lines.join(" · ") || tr("unknown");
}

function renderFamilyQueue() {
  const field = $("familyQueueField");
  const select = $("familyQueue");
  field.hidden = !state.familyQa.enabled;
  if (!state.familyQa.enabled) return;
  select.innerHTML = state.familyQa.families.map((family) =>
    `<option value="${escapeHtml(family.family_id)}">${escapeHtml(familyOptionLabel(family))}</option>`,
  ).join("");
  select.value = state.activeFamilyId;
}

function currentRepresentative() {
  const family = activeFamily();
  return family?.representatives.find((rep) => rep.case_id === state.current?.case_id) || null;
}

function renderFamilyQaCard() {
  const card = $("familyQaCard");
  const family = activeFamily();
  const rep = currentRepresentative();
  card.hidden = !family;
  document.body.classList.toggle("family-qa-mode", Boolean(family));
  if (!family) return;
  const index = Math.max(0, family.representatives.findIndex((item) => item.case_id === rep?.case_id));
  $("familyQaTitle").textContent = `${family.family_id} · ${family.family_size} cases`;
  $("familyQaRep").textContent = `rep ${index + 1} / ${family.representatives.length}`;
  const progress = state.familyQa.progress || {};
  $("familyQaProgress").textContent = `${progress.reviewed_families || 0} / ${progress.total_families || 0} families reviewed · ${progress.approved_cases || 0} / ${progress.total_cases || 0} cases approved`;
  $("familyCalibration").textContent = rep?.calibration_relation || family.calibration_relation || "—";
  $("familyMetalState").textContent = `${rep?.candidate_metal_state || "—"} → ${rep?.reference_metal_state || "—"}`;
  $("familyRepeat").textContent = rep?.repeat_count || "—";
  $("familyMapping").textContent = rep?.mapping_source || "—";
  $("familyTransformation").textContent = transformationSummary(
    rep?.canonical_transformation || family.canonical_transformation,
    rep?.repeat_count || family.repeat_count,
  );
  document.querySelectorAll("[data-rep-mark]").forEach((button) => {
    button.classList.toggle("selected", button.dataset.repMark === rep?.qa_mark);
  });
  document.querySelectorAll("[data-family-decision]").forEach((button) => {
    button.classList.toggle("selected", button.dataset.familyDecision === family.decision);
  });
}

async function loadFamilyQa({ preserveFamily = true } = {}) {
  const data = await api("/api/family-qa");
  state.familyQa = data;
  if (data.enabled && (!preserveFamily || !data.families.some((family) => family.family_id === state.activeFamilyId))) {
    state.activeFamilyId = data.families.find((family) => !family.decision)?.family_id || data.families[0]?.family_id || "";
  }
  renderFamilyQueue();
  renderFamilyQaCard();
}

async function mutateFamilyQa(action, value) {
  const family = activeFamily();
  const rep = currentRepresentative();
  if (!family) return;
  const result = await api("/api/family-qa", {
    method: "POST",
    body: JSON.stringify({ family_id: family.family_id, case_id: rep?.case_id || "", action, value }),
  });
  if (result.mutation) {
    state.reviewHistory.unshift(result.mutation);
    state.reviewHistory = state.reviewHistory.slice(0, 20);
  }
  state.familyQa = result.family_qa;
  renderFamilyQueue();
  renderFamilyQaCard();
  renderReviewHistory();
}

async function loadCases(reset = false) {
  const family = activeFamily();
  if (family) {
    state.offset = 0;
    state.cases = family.representatives.map((rep) => ({ ...rep.case, family_rep: rep }));
    state.total = state.cases.length;
    renderCaseList();
    return;
  }
  if (reset) state.offset = 0;
  const params = new URLSearchParams();
  const category = $("categoryFilter").value;
  const status = $("statusFilter").value;
  const reviewReason = $("reviewReasonFilter").value;
  const q = $("searchBox").value.trim();
  const triageBucket = $("triageFilter").value;
  if (category) params.set("category", category);
  if (status) params.set("status", status);
  if (reviewReason) params.set("reviewer", reviewReason);
  if (q) params.set("q", q);
  if (triageBucket) params.set("triage_bucket", triageBucket);
  params.set("limit", state.limit);
  params.set("offset", state.offset);
  const data = await api(`/api/cases?${params.toString()}`);
  state.cases = data.items || [];
  state.total = data.total || 0;
  renderTriageBucketOptions(data.triage_bucket_counts || {});
  renderCaseList();
}

function renderCaseList() {
  const filteredTriageBucket = $("triageFilter").value;
  $("caseList").innerHTML = state.cases
    .map((item) => {
      const selected = state.current && state.current.case_id === item.case_id ? "selected" : "";
      const status = statusLabel(item.review_status);
      const category = hasValue(item.category)
        ? badge(categoryLabel(item.category), `queue-tag ${categoryKind(item.category)}`)
        : "";
      const fixture = item.fixture
        ? badge(msg("fixture", { kind: item.fixture.kind }), "fixture-tag")
        : "";
      const triage =
        item.triage_bucket && item.triage_bucket !== filteredTriageBucket
          ? badge(triageBucketLabel(item.triage_bucket), "triage-tag")
          : "";
      return `
        <button class="case-item ${selected}" data-case-id="${escapeHtml(item.case_id)}" type="button">
          <span class="row"><strong>${escapeHtml(item.case_id)}</strong><span>#${item.row_index}</span></span>
          <span class="row">
            ${category}
            ${badge(status, `review-tag ${item.review_status ? reviewStatusKind(item.review_status) : ""}`)}
            ${triage}
            ${fixture}
          </span>
        </button>`;
    })
    .join("");
  $("pageInfo").textContent = `${state.offset + 1}-${Math.min(state.offset + state.limit, state.total)} / ${state.total}`;
  document.querySelectorAll(".case-item").forEach((button) => {
    button.addEventListener("click", () => loadCase(button.dataset.caseId));
  });
}

function isCurrentCaseRequest(token, caseId) {
  return token === state.caseRequestToken && state.current?.case_id === caseId;
}

async function loadCase(caseId) {
  const token = ++state.caseRequestToken;
  let item;
  try {
    item = await api(`/api/cases/${encodeURIComponent(caseId)}`);
  } catch (error) {
    if (token !== state.caseRequestToken) return;
    throw error;
  }
  if (token !== state.caseRequestToken) return;
  state.current = item;
  state.twoDMode = hasHydrogenAssignment(item) ? "hydrogen" : "skeleton";
  syncTwoDModeButtons();
  const url = new URL(window.location.href);
  url.searchParams.set("case", item.case_id);
  window.history.replaceState({}, "", url);
  state.currentLiveCandidate = null;
  state.xyzViewer = null;
  state.referenceXyzViewer = null;
  state.referenceXyzSdf = "";
  state.referenceXyzFailure = "";
  state.mappedComparison = {};
  renderMappedComparisonNote();
  state.referenceRenderStatus = "render_failed";
  state.referenceRenderError = "";
  state.graphEvidence = {};
  state.graphEvidenceLoadingCaseId = "";
  state.xyzLoadStatus = "unknown";
  state.xyzLoadError = "";
  renderCaseHeader();
  renderReviewerSummary();
  populateReviewForm();
  renderVersionComparison();
  renderDiagnostics();
  renderReviewerDetails();
  renderFamilyQaCard();
  await Promise.all([
    loadXyz(item, token),
    loadReferenceXyz(item, token),
    loadCandidateSdf(item, token),
  ]);
  if (!isCurrentCaseRequest(token, item.case_id)) return;
  await loadPair(token);
  if (!isCurrentCaseRequest(token, item.case_id)) return;
  if ($("reviewerDetails").open) await loadGraphEvidence(token);
  if (!isCurrentCaseRequest(token, item.case_id)) return;
  renderCaseList();
}

function hasHydrogenAssignment(item) {
  const triage = item?.triage;
  if (!triage) return false;
  return parseJsonArray(triage.hydrogen_assignment_diff).length > 0
    || String(triage.reason_tags || "").split(/[;,|]/).some((value) =>
      ["h_assignment", "hydrogen_assignment", "hydrogen-assignment"].includes(value.trim()),
    );
}

function syncTwoDModeButtons() {
  document.querySelectorAll(".two-d-mode").forEach((button) => {
    button.classList.toggle("active", button.dataset.mode === state.twoDMode);
  });
}

function renderCaseHeader() {
  const item = state.current;
  const openTrace = $("openTrace");
  openTrace.disabled = !item;
  openTrace.title = item ? msg("traceTitle", { caseId: item.case_id }) : tr("selectCaseFirst");
  const fixture = item.fixture;
  $("caseTitle").textContent = item.case_id;
  const context = hasValue(item.row_index) ? [`${tr("row")} ${item.row_index}`] : [];
  $("caseContext").textContent = context.length ? context.join(" · ") : tr("notProvided");
  const candidateStatus = item.live_candidate_status || item.candidate_status;
  const headerItems = [];
  if (hasValue(item.category)) {
    headerItems.push(badge(categoryLabel(item.category), `queue-tag ${categoryKind(item.category)}`));
  }
  if (hasValue(item.triage_bucket)) {
    headerItems.push(badge(triageBucketLabel(item.triage_bucket), "triage-tag"));
  }
  headerItems.push(
    item.review_status
      ? badge(statusLabel(item.review_status), `review-tag ${reviewStatusKind(item.review_status)}`)
      : badge(tr("unreviewed"), "review-tag"),
  );
  if (fixture) headerItems.push(badge(msg("fixture", { kind: fixture.kind }), "fixture-tag"));
  $("caseMeta").innerHTML = headerItems.join(" ");
  const audit = [];
  if (hasValue(item.source)) audit.push(item.source);
  if (candidateStatus && candidateStatus !== "ok") audit.push(`${tr("currentCandidateLabel")}: ${candidateStatus}`);
  if (item.review_status && hasValue(item.reviewer)) {
    audit.push(msg("reviewedBy", { reviewer: item.reviewer }));
  }
  if (item.review_status && hasValue(item.updated_at)) {
    audit.push(msg("updatedAt", { updatedAt: formatTimestamp(item.updated_at) }));
  }
  $("reviewAudit").textContent = audit.join(" · ");
  const removeFixture = $("removeFixture");
  removeFixture.hidden = !fixture;
  removeFixture.disabled = !fixture;
  removeFixture.title = fixture
    ? msg("removeFixtureTitle", { file: fixture.structure_file, caseId: item.case_id })
    : tr("noFixture");
}

function hasValue(value) {
  return value !== null && value !== undefined && String(value).trim() !== "";
}

function formatTimestamp(value) {
  if (!hasValue(value)) return tr("notProvided");
  const date = new Date(value);
  return Number.isNaN(date.getTime()) ? String(value) : date.toLocaleString();
}

function referenceState(item) {
  if (!hasValue(item?.reference_smiles)) return "missing";
  const parseStatus = String(item.reference_parse_status || "").toLowerCase();
  if (parseStatus && parseStatus !== "ok") return "parse_invalid";
  return state.referenceRenderStatus === "available" ? "available" : "render_failed";
}

function referenceLabel(item) {
  const status = referenceState(item);
  if (status === "missing") return tr("missing");
  if (status === "parse_invalid") return tr("invalid");
  if (status === "render_failed") return tr("referenceRenderFailed");
  return tr("valid");
}

function summaryValue(value) {
  return hasValue(value) ? String(value) : tr("notProvided");
}

function formulaSummary(item) {
  const fields = [
    formulaDisplayStatus(item),
    item.reference_formula_match,
    item.xyz_formula,
    item.reference_formula_with_h,
  ].filter(hasValue);
  return fields.length ? fields.join(" · ") : tr("notProvided");
}

function formulaDisplayStatus(item) {
  const status = String(item.reference_formula_check_status || "").toLowerCase();
  if (status === "comparison_skipped" && ["equivalence_timeout", "candidate_reparse_failure"]
    .includes(item.reference_diagnostic_group)) return "ok";
  return item.reference_formula_check_status;
}

function formulaCompact(item) {
  const match = item.reference_formula_match;
  const status = String(item.reference_formula_check_status || "").toLowerCase();
  if (!hasValue(match) && !status) return tr("statusUnknownCompact");
  if (status === "comparison_skipped") return tr("statusUnknownCompact");
  if (match === "True" || status === "ok") return tr("formulaOk");
  if (["not_applicable", "missing_reference_smiles"].includes(status)) return tr("statusUnknownCompact");
  return summaryValue(item.reference_formula_check_status || match);
}

function assessmentSummary(item) {
  const fields = [
    item.accuracy_assessment_status,
    item.tmqmg_answer_assessment,
    item.molgr_answer_assessment,
  ].filter(hasValue);
  return fields.length ? fields.join(" · ") : tr("notProvided");
}

function snapshotSummary(item) {
  if (item.live_matches_candidate_snapshot === true) return tr("generatedSnapshotMatch");
  if (item.live_matches_candidate_snapshot === false) return tr("generatedSnapshotMismatch");
  if (hasValue(item.live_candidate_status)) return tr("generatedSnapshotIncomparable");
  return tr("notProvided");
}

function snapshotCompact(item) {
  if (item.live_matches_candidate_snapshot === true) return tr("snapshotCurrent");
  if (item.live_matches_candidate_snapshot === false) return tr("snapshotDifferent");
  return tr("snapshotUnknown");
}

function referenceCompact(item) {
  const status = referenceState(item);
  if (status === "available") return tr("statusOk");
  if (status === "missing") return tr("statusMissingCompact");
  if (status === "parse_invalid") return tr("statusInvalidCompact");
  return tr("renderFailedCompact");
}

function candidateCompact(item) {
  const status = String(item.live_candidate_status || item.candidate_status || "").toLowerCase();
  if (status === "ok") return tr("statusOk");
  if (["failed", "error", "unavailable"].includes(status)) return tr("statusUnavailableCompact");
  return status || tr("statusUnknownCompact");
}

function assessabilityCompact(item) {
  const value = item.accuracy_assessment_status || item.tmqmg_answer_assessment || item.molgr_answer_assessment;
  if (!hasValue(value)) return tr("statusUnknownCompact");
  return String(value).toLowerCase() === "assessable" ? tr("statusOk") : String(value);
}

function reviewFocusKey(item) {
  if (["failed", "error", "unavailable"].includes(String(item.live_candidate_status || item.candidate_status || "").toLowerCase())) {
    return "focusCandidateFailure";
  }
  if (["missing", "parse_invalid", "render_failed"].includes(referenceState(item))) return "focusReference";
  if (["equivalence_timeout", "candidate_reparse_failure", "comparison_skipped"]
    .includes(item.reference_diagnostic_group)) return "focusReferenceComparison";
  const formulaStatus = String(item.reference_formula_check_status || "").toLowerCase();
  if ((item.reference_formula_match === "False" && formulaStatus !== "comparison_skipped")
    || (formulaStatus && !["ok", "not_applicable", "comparison_skipped"].includes(formulaStatus))) {
    return "focusFormula";
  }
  const assessability = assessmentSummary(item).toLowerCase();
  if (assessability !== tr("notProvided").toLowerCase() && /not_assessable|unassessable|limited/.test(assessability)) {
    return "focusAssessment";
  }
  if (item.live_matches_candidate_snapshot === false) return "focusSnapshot";
  return "focusComparison";
}

function renderReviewerSummary() {
  const item = state.current;
  if (!item) return;
  const referenceStatus = referenceState(item);
  const candidateStatus = String(item.live_candidate_status || item.candidate_status || "").toLowerCase();
  const snapshotDifferent = item.live_matches_candidate_snapshot === false;
  const snapshotKnown = typeof item.live_matches_candidate_snapshot === "boolean";
  const provenance = snapshotDifferent ? compactProvenanceReason(item.live_candidate_equivalence_reason) : "";
  const secondary = [
    ["q", summaryValue(item.total_charge), "qTooltip"],
    ["M", summaryValue(item.spin_multiplicity), "multiplicityTooltip"],
    ["radicals", summaryValue(item.total_radical_electrons), "radicalsTooltip"],
    [tr("formulaLabel"), formulaCompact(item), "formulaTooltip"],
  ];
  const snapshotText = snapshotKnown
    ? snapshotDifferent
      ? `<span class="status-pill drift">${escapeHtml(tr("differentCompact"))}</span>`
      : `<span class="snapshot-same">${escapeHtml(tr("sameCompact"))}</span>`
    : escapeHtml(tr("statusUnknownCompact"));
  const mainStatus = [
    `${escapeHtml(tr("currentCandidateLabel"))} ${candidateStatus === "ok" ? "✓" : escapeHtml(candidateStatus || tr("statusUnknownCompact"))}`,
    `${escapeHtml(tr("summaryReference"))} ${referenceStatus === "available" ? "✓" : escapeHtml(referenceStatus === "missing" ? tr("missingCompact") : referenceStatus === "parse_invalid" ? tr("invalidCompact") : tr("renderFailedCompact"))}`,
    `${escapeHtml(tr("currentVsSnapshotLabel"))} ${snapshotText}`,
  ];
  $("reviewSummary").innerHTML = `
    <div class="status-line">${mainStatus.join('<span class="status-separator" aria-hidden="true"> · </span>')}</div>
    ${provenance ? `<div class="drift-reason">${escapeHtml(provenance)}</div>` : ""}
    <div class="status-line metadata-line">${secondary
      .map(
        ([label, value, tooltip]) =>
          `<span title="${escapeHtml(tr(tooltip))}">${escapeHtml(label)} ${escapeHtml(value)}</span>`,
      )
      .join('<span class="status-separator" aria-hidden="true"> · </span>')}</div>`;
  renderTriageEvidence(item);
  updateReferenceVisual();
}

function parseJsonArray(value) {
  if (Array.isArray(value)) return value;
  if (!hasValue(value)) return [];
  try {
    const parsed = JSON.parse(value);
    return Array.isArray(parsed) ? parsed : [];
  } catch (error) {
    return [];
  }
}

const TRIAGE_BUCKET_LABELS = {
  zh: {
    strong_xyz_candidate_evidence: "XYZ→候选",
    strong_xyz_reference_evidence: "XYZ→参考",
    possible_redox_representation: "氧化态/表示",
    metal_coordination_ambiguous: "配位待判",
    organic_topology_manual: "有机拓扑",
    complex_multi_difference: "多重差异",
    reference_integrity_issue: "参考问题",
    unknown: "未分类",
  },
  en: {
    strong_xyz_candidate_evidence: "XYZ→Candidate",
    strong_xyz_reference_evidence: "XYZ→Reference",
    possible_redox_representation: "Oxidation/representation",
    metal_coordination_ambiguous: "Coordination review",
    organic_topology_manual: "Organic topology",
    complex_multi_difference: "Multiple differences",
    reference_integrity_issue: "Reference issue",
    unknown: "Unclassified",
  },
};

function triageBucketLabel(bucket) {
  return TRIAGE_BUCKET_LABELS[state.language]?.[bucket] || bucket;
}

function localizeTriageFilterOptions() {
  document.querySelectorAll("#triageFilter option[value]").forEach((option) => {
    if (!option.value) {
      option.textContent = tr("all");
      return;
    }
    const count = option.dataset.count;
    option.textContent = `${triageBucketLabel(option.value)}${count ? ` (${count})` : ""}`;
  });
}

function mappingLabel(confidence) {
  if (confidence === "unique_graph_mapping") return "unique";
  return confidence || "—";
}

function evidenceRow(label, value, mono = false) {
  if (!hasValue(value)) return "";
  const valueClass = mono ? ' class="mono"' : "";
  return `<div class="triage-evidence-row"><span>${escapeHtml(label)}</span><strong${valueClass}>${escapeHtml(value)}</strong></div>`;
}

function metalEdgeEvidence(triage, edge) {
  const elements = Array.isArray(edge.elements) ? edge.elements : [];
  const atoms = Array.isArray(edge.candidate_atoms) ? edge.candidate_atoms : [];
  const metalSymbols = String(triage.metal_elements || "")
    .split(/[|,;\s]+/)
    .filter(Boolean);
  let metalIndex = elements.findIndex((element) => metalSymbols.includes(String(element)));
  if (metalIndex < 0) metalIndex = 0;
  const ligandIndex = metalIndex === 0 ? 1 : 0;
  const side = edge.edge_present_in;
  const difference = side === "reference"
    ? "Reference-only coordination"
    : side === "candidate"
      ? "Candidate-only coordination"
      : "Coordination difference";
  const distance = Number(edge.distance);
  const candidatePresence =
    side === "candidate" ? "present" : side === "reference" ? "absent" : "—";
  const referencePresence =
    side === "reference" ? "present" : side === "candidate" ? "absent" : "—";
  const cn =
    side === "reference" ? edge.reference_coordination_number : edge.candidate_coordination_number;
  const shell = edge.inside_agreed_shell_range === true ? "disputed atom in donor shell" : "";
  return [
    evidenceRow("差异", difference),
    evidenceRow(
      "原子",
      `${elements[metalIndex] || "?"} · XYZ #${atoms[metalIndex] ?? "?"} ↔ ${elements[ligandIndex] || "?"} · XYZ #${atoms[ligandIndex] ?? "?"}`,
      true,
    ),
    evidenceRow("XYZ 距离", Number.isFinite(distance) ? `${distance.toFixed(3)} Å` : "—", true),
    evidenceRow("Candidate", candidatePresence),
    evidenceRow("Reference", referencePresence),
    evidenceRow("Mapping", mappingLabel(triage.mapping_confidence), true),
    evidenceRow("Trace", [hasValue(cn) ? `CN=${cn}` : "", shell].filter(Boolean).join(" · ")),
  ].join("");
}

function hydrogenEvidence(triage, hydrogen) {
  const candidateDistance = Number(hydrogen.candidate_distance);
  const referenceDistance = Number(hydrogen.reference_distance);
  const margin = Number(hydrogen.distance_margin);
  const distanceParts = [];
  if (Number.isFinite(candidateDistance)) {
    distanceParts.push(
      `${hydrogen.candidate_center_element || "?"}-H ${candidateDistance.toFixed(3)} Å`,
    );
  }
  if (Number.isFinite(referenceDistance)) {
    distanceParts.push(
      `${hydrogen.reference_center_element || "?"}-H ${referenceDistance.toFixed(3)} Å`,
    );
  }
  return [
    evidenceRow("差异", "H assignment"),
    evidenceRow("H", `H · XYZ #${hydrogen.h_atom ?? "?"}`, true),
    evidenceRow(
      "Candidate",
      `${hydrogen.candidate_center_element || "?"} #${hydrogen.candidate_center ?? "?"}`,
      true,
    ),
    evidenceRow(
      "Reference",
      `${hydrogen.reference_center_element || "?"} #${hydrogen.reference_center ?? "?"}`,
      true,
    ),
    evidenceRow("XYZ 距离", distanceParts.join(" · "), true),
    evidenceRow("Margin", Number.isFinite(margin) ? `${margin.toFixed(3)} Å` : "—", true),
    evidenceRow("Mapping", mappingLabel(triage.mapping_confidence), true),
  ].join("");
}

function redoxEvidence(triage) {
  const coordination = parseJsonArray(triage.metal_coordination_diff);
  const coordinationText = coordination.length
    ? `${coordination.length} disputed edge${coordination.length === 1 ? "" : "s"}`
    : "none";
  return [
    evidenceRow("差异", "Metal oxidation state / representation"),
    evidenceRow("Metal", triage.metal_elements, true),
    evidenceRow("Candidate metal", triage.candidate_metal_state, true),
    evidenceRow("Reference metal", triage.reference_metal_state, true),
    evidenceRow("Metal Δ", triage.metal_charge_delta, true),
    evidenceRow(
      "Ligand compensation",
      hasValue(triage.ligand_charge_delta)
        ? `Δ ${triage.ligand_charge_delta} · Candidate ${triage.candidate_ligand_charge_sum || "—"} · Reference ${triage.reference_ligand_charge_sum || "—"}`
        : "—",
      true,
    ),
    evidenceRow("Coordination Δ", coordinationText),
    evidenceRow("Mapping", mappingLabel(triage.mapping_confidence), true),
  ].join("");
}

function referenceProblemEvidence(item) {
  const group = item.reference_diagnostic_group || item.triage?.reference_diagnostic_group;
  const reason = item.reference_diagnostic_reason || item.triage?.reference_diagnostic_reason;
  return [
    evidenceRow(tr("diagnosticLabel"), tr(`referenceDiagnostics.${group}`, group)),
    evidenceRow(tr("diagnosticReasonLabel"), reason),
  ].join("");
}

function renderTriageEvidence(item) {
  const panel = $("triageEvidence");
  const triage = item?.triage;
  if (!triage) {
    panel.hidden = true;
    panel.innerHTML = "";
    return;
  }
  const metalEdges = parseJsonArray(triage.metal_coordination_diff);
  const hydrogenAssignments = parseJsonArray(triage.hydrogen_assignment_diff);
  let evidence = "";
  if (triage.triage_bucket === "reference_integrity_issue") {
    evidence = referenceProblemEvidence(item);
  } else if (triage.triage_bucket === "possible_redox_representation") {
    evidence = redoxEvidence(triage);
  } else if (hydrogenAssignments.length) {
    evidence = hydrogenAssignments
      .map((entry) => hydrogenEvidence(triage, entry))
      .join('<div class="triage-evidence-separator"></div>');
    if (metalEdges.length) {
      evidence += metalEdges.map((edge) => metalEdgeEvidence(triage, edge)).join("");
    }
  } else if (metalEdges.length) {
    evidence = metalEdges
      .map((edge) => metalEdgeEvidence(triage, edge))
      .join('<div class="triage-evidence-separator"></div>');
  } else {
    evidence = evidenceRow("Evidence", triage.xyz_evidence_summary || triage.machine_reason || "—");
  }
  panel.hidden = false;
  panel.innerHTML = `
    <div class="triage-evidence-head"><strong>${escapeHtml(tr("triageEvidence"))}</strong>
      <span>${escapeHtml(triageBucketLabel(triage.triage_bucket || ""))}</span>
    </div>
    <div class="triage-evidence-grid">${evidence}</div>
    ${hasValue(triage.trace_evidence_summary) ? `<details class="triage-trace"><summary>${escapeHtml(tr("fullTraceEvidence"))}</summary><div>${escapeHtml(triage.trace_evidence_summary)}</div></details>` : ""}`;
}

function compactProvenanceReason(reason) {
  return String(reason || "")
    .replace(/^Not equivalent:\s*/i, "")
    .replace(/\.$/, "");
}

function updateReferenceVisual() {
  const item = state.current;
  if (!item) return;
  const status = referenceState(item);
  $("primaryVisual").hidden = false;
  $("secondaryVisual").hidden = false;
  if (status === "missing") showReferenceMessage(tr("referenceMissingNotice"));
  if (status === "parse_invalid") showReferenceMessage(tr("referenceInvalidNotice"));
  if (status === "render_failed") {
    showReferenceMessage(tr("referenceRenderFailed"), state.referenceRenderError);
  }
}

function populateReviewForm() {
  const item = state.current;
  $("correctedSmiles").value = item.corrected_smiles || "";
  $("correctedMolblock").value = item.corrected_molblock || "";
  $("notes").value = item.notes || "";
  $("reviewer").value = item.reviewer || "";
  document.querySelectorAll(".decision").forEach((button) => {
    button.classList.toggle("selected", button.dataset.status === item.review_status);
  });
  if (item.review_status === "manual_reference") {
    $("manualEditorDetails").open = true;
  }
}

function renderReviewerDetails() {
  const item = state.current;
  if (!item) return;
  const smilesRows = [
    [tr("currentCandidateSmiles"), item.live_candidate_smiles],
    [tr("candidateSnapshotSmiles"), item.candidate_snapshot_smiles],
    [tr("referenceSmiles"), item.reference_smiles],
  ];
  $("reviewerSmiles").innerHTML = smilesRows
    .map(
      ([label, value]) =>
        `<dt>${escapeHtml(label)}</dt><dd><code>${escapeHtml(hasValue(value) ? value : tr("notProvided"))}</code></dd>`,
    )
    .join("");
  renderDiagnosticList(
    "reviewerProvenance",
    [
      [tr("currentCandidateStatus"), item.live_candidate_status || item.candidate_status],
      [tr("currentVsSnapshot"), snapshotSummary(item)],
      [tr("equivalenceMethod"), item.live_candidate_equivalence_method],
      [tr("equivalenceReason"), item.live_candidate_equivalence_reason],
      [tr("snapshotRuntime"), item.candidate_snapshot_runtime],
    ],
    false,
  );
  renderGraphEvidence();
}

function jumpToReviewerDetails() {
  const details = $("reviewerDetails");
  details.open = true;
  loadGraphEvidence(state.caseRequestToken);
  requestAnimationFrame(() => details.scrollIntoView({ behavior: "smooth", block: "start" }));
}

function signedCharge(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return tr("notProvided");
  return number > 0 ? `+${number}` : String(number);
}

function metalSummary(metals) {
  if (!Array.isArray(metals) || !metals.length) return "—";
  return metals
    .map((metal) => `${metal.element}${metal.index} ${signedCharge(metal.formal_charge)}`)
    .join(", ");
}

function renderGraphEvidence() {
  const metrics = [
    ["totalFormalCharge", "total_formal_charge", (value) => signedCharge(value)],
    ["atomCount", "atom_count", String],
    ["explicitHCount", "explicit_h_count", String],
    ["totalRadicalElectrons", "total_radical_electrons", String],
    ["metals", "metals", metalSummary],
  ];
  const candidate = state.graphEvidence.candidate?.summary;
  const reference = state.graphEvidence.reference?.summary;
  $("graphSummaryRows").innerHTML = metrics
    .map(([labelKey, field, format]) => {
      const candidateValue = candidate && candidate[field] !== undefined
        ? format(candidate[field])
        : "—";
      const referenceValue = reference && reference[field] !== undefined
        ? format(reference[field])
        : "—";
      return `<tr><th>${escapeHtml(tr(labelKey))}</th><td>${escapeHtml(candidateValue)}</td><td>${escapeHtml(referenceValue)}</td></tr>`;
    })
    .join("");
  renderGraphInspector("candidate");
  renderGraphInspector("reference");
}

function relevantAtomIndices(graph, showFull) {
  if (!graph || !Array.isArray(graph.atoms)) return new Set();
  if (showFull) return new Set(graph.atoms.map((atom) => atom.index));
  const metalIndices = new Set(
    graph.atoms.filter((atom) => atom.is_metal).map((atom) => atom.index),
  );
  return new Set(
    graph.atoms
      .filter(
        (atom) =>
          !["C", "H"].includes(atom.element) ||
          atom.is_metal ||
          atom.neighbours.some((neighbour) => metalIndices.has(neighbour.index)),
      )
      .map((atom) => atom.index),
  );
}

function atomLabel(element, index) {
  return `${element}${index}`;
}

function renderGraphInspector(kind) {
  const graph = state.graphEvidence[kind];
  const prefix = kind === "candidate" ? "candidate" : "reference";
  const atomRows = $(`${prefix}AtomRows`);
  const bondRows = $(`${prefix}BondRows`);
  const error = $(`${prefix}GraphError`);
  if (!graph) {
    atomRows.innerHTML = `<tr><td colspan="7">${escapeHtml(tr("graphLoading"))}</td></tr>`;
    bondRows.innerHTML = `<tr><td colspan="3">${escapeHtml(tr("graphLoading"))}</td></tr>`;
    error.hidden = true;
    error.textContent = "";
    return;
  }
  if (graph.error) {
    atomRows.innerHTML = `<tr><td colspan="7">${escapeHtml(tr("graphUnavailable"))}</td></tr>`;
    bondRows.innerHTML = `<tr><td colspan="3">${escapeHtml(tr("graphUnavailable"))}</td></tr>`;
    error.hidden = false;
    error.textContent = tr("graphUnavailable");
    return;
  }
  error.hidden = true;
  error.textContent = "";
  const showFull = document.querySelector(`.show-full-graph[data-kind="${kind}"]`)?.checked;
  const visible = relevantAtomIndices(graph, showFull);
  const atoms = graph.atoms.filter((atom) => visible.has(atom.index));
  atomRows.innerHTML = atoms.length
    ? atoms
        .map(
          (atom) => `<tr>
            <td>${atom.index}</td><td>${escapeHtml(atom.element)}</td>
            <td>${escapeHtml(signedCharge(atom.formal_charge))}</td>
            <td>${atom.radical_electrons}</td><td>${atom.explicit_h ?? "—"}</td>
            <td>${atom.implicit_h ?? "—"}</td>
            <td>${escapeHtml(atom.neighbours.map((neighbour) => atomLabel(neighbour.element, neighbour.index)).join(", "))}</td>
          </tr>`,
        )
        .join("")
    : `<tr><td colspan="7">${escapeHtml(tr("noRelevantAtoms"))}</td></tr>`;
  const bonds = graph.bonds.filter(
    (bond) => showFull || visible.has(bond.begin_atom) || visible.has(bond.end_atom),
  );
  bondRows.innerHTML = bonds.length
    ? bonds
        .map((bond) => {
          const begin = atomLabel(bond.begin_element, bond.begin_atom);
          const end = atomLabel(bond.end_element, bond.end_atom);
          const connector = bond.directional
            ? "→"
            : bond.type === "double"
              ? "="
              : bond.type === "triple"
                ? "≡"
                : bond.type === "aromatic"
                  ? "↔"
                  : "–";
          return `<tr><td>${bond.index}</td><td>${escapeHtml(`${begin} ${connector} ${end}`)}</td><td>${escapeHtml(bond.type)}</td></tr>`;
        })
        .join("")
    : `<tr><td colspan="3">—</td></tr>`;
}

async function loadGraphEvidence(token = state.caseRequestToken) {
  const item = state.current;
  if (!item || !isCurrentCaseRequest(token, item.case_id)) return;
  if (state.graphEvidenceLoadingCaseId === item.case_id) return;
  if (state.graphEvidence.candidate && state.graphEvidence.reference) return;
  state.graphEvidenceLoadingCaseId = item.case_id;
  renderGraphEvidence();
  const kinds = ["candidate", "reference"];
  const results = await Promise.all(
    kinds.map(async (kind) => {
      try {
        return await api(
          `/api/cases/${encodeURIComponent(item.case_id)}/graph?kind=${encodeURIComponent(kind)}`,
        );
      } catch (error) {
        return { kind, error: error.message };
      }
    }),
  );
  if (!isCurrentCaseRequest(token, item.case_id)) return;
  state.graphEvidence = Object.fromEntries(results.map((result) => [result.kind, result]));
  state.graphEvidenceLoadingCaseId = "";
  renderGraphEvidence();
}

function renderDiagnostics() {
  const item = state.current;
  const smilesPairs = [
    ["reference_smiles", item.reference_smiles],
    ["candidate_smiles", item.candidate_smiles],
    ["candidate_snapshot_smiles", item.candidate_snapshot_smiles],
    ["live_candidate_smiles", item.live_candidate_smiles],
    ["candidate_organic", item.candidate_organic_smiles],
    ["reference_organic", item.reference_organic_smiles],
  ];
  const assessmentPairs = [
    ["reference_formula_status", formulaDisplayStatus(item)],
    ["reference_diagnostic_group", item.reference_diagnostic_group],
    ["reference_diagnostic_reason", item.reference_diagnostic_reason],
    ["xyz_formula", item.xyz_formula],
    ["reference_formula_with_h", item.reference_formula_with_h],
    ["reference_formula_mismatch", item.reference_formula_mismatch_detail],
    ["reference_answer_status", item.reference_answer_status],
    ["reference_answer_reason", item.reference_answer_reason],
    ["accuracy_assessment_status", item.accuracy_assessment_status],
    ["accuracy_assessment_reason", item.accuracy_assessment_reason],
    ["tmqmg_answer_assessment", item.tmqmg_answer_assessment],
    ["molgr_answer_assessment", item.molgr_answer_assessment],
    ["error", item.error],
  ];
  const runtimePairs = [
    ["candidate_snapshot_runtime", item.candidate_snapshot_runtime],
    ["live_candidate_smiles_exact_match", item.live_candidate_smiles_exact_match],
    ["live_matches_candidate_snapshot", item.live_matches_candidate_snapshot],
    ["live_candidate_reason", item.live_candidate_equivalence_reason],
    ["reference_xyz_failure", state.referenceXyzFailure],
  ];
  renderDiagnosticList("assessmentDiagnostics", assessmentPairs);
  renderDiagnosticList("runtimeDiagnostics", runtimePairs);

  const represented = new Set([
    "fixture", "case_id", "row_index", "source", "category", "xyz_path",
    "total_charge", "total_radical_electrons", "spin_multiplicity", "reference_smiles",
    "candidate_smiles", "candidate_organic_smiles", "reference_organic_smiles",
    "candidate_status", "review_status", "corrected_smiles", "corrected_molblock", "notes",
    "reviewer", "updated_at", "candidate_snapshot_runtime", "candidate_snapshot_smiles",
    "candidate_snapshot_status", "live_candidate_status", "live_candidate_smiles",
    "live_candidate_smiles_exact_match", "live_matches_candidate_snapshot",
    "live_candidate_equivalence_method", "live_candidate_equivalence_reason",
    "reference_formula_check_status", "xyz_formula", "reference_formula_with_h",
    "reference_diagnostic_group", "reference_diagnostic_reason",
    "reference_formula_mismatch_detail", "reference_answer_status", "reference_answer_reason",
    "accuracy_assessment_status", "accuracy_assessment_reason", "tmqmg_answer_assessment",
    "molgr_answer_assessment", "error",
  ]);
  const metadataPairs = Object.entries(item)
    .filter(([key, value]) => !represented.has(key) && hasValue(value))
    .sort(([left], [right]) => left.localeCompare(right));
  renderDiagnosticList("metadataDiagnostics", metadataPairs, false);

  $("diagnostics").innerHTML = [...smilesPairs, ...assessmentPairs, ...runtimePairs]
    .map(([key, value]) => {
      const displayValue = value === null || value === undefined ? "" : String(value);
      return `<dt>${escapeHtml(tr(`diagnostics.${key}`, key))}</dt><dd>${escapeHtml(displayValue)}</dd>`;
    })
    .join("");
}

function renderDiagnosticList(id, pairs, translateKeys = true) {
  const available = pairs.filter(([, value]) => hasValue(value));
  $(id).innerHTML = available.length
    ? available
        .map(([key, value]) => {
          const label = translateKeys ? tr(`diagnostics.${key}`, key) : key;
          const display = typeof value === "object" ? JSON.stringify(value, null, 2) : String(value);
          return `<dt>${escapeHtml(label)}</dt><dd>${escapeHtml(display)}</dd>`;
        })
        .join("")
    : `<div class="details-empty">${escapeHtml(tr("notProvided"))}</div>`;
}

function statusBadge(status) {
  if (!status) return badge("missing", "bad");
  if (status === "ok") return badge(status, "ok");
  if (status === "not_run") return badge(status);
  return badge(status, "bad");
}

function versionResultCard(label, status) {
  return `
    <article class="version-card">
      <header>
        <strong>${escapeHtml(label)}</strong>
        ${statusBadge(status)}
      </header>
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
      <h4>${escapeHtml(tr("versionComparison"))}</h4>
      ${badge(mismatch ? tr("disagreement") : tr("consistent"), mismatch ? "warn" : "ok")}
    </header>
    <div class="version-grid">
      ${versionResultCard("Python 3.8 · candidate_cpp", py38Status)}
      ${versionResultCard("Python 3.10 · candidate_cpp", py310Status)}
    </div>`;
}

function renderCandidateSdfStatus() {
  const item = state.current;
  if (!item) return;
  const status = $("candidateSdfStatus");
  const modelStatus = $("candidateModelStatus");
  if (item.live_candidate_status === "error") {
    status.textContent = tr("error");
    modelStatus.textContent = tr("error");
    return;
  }
  if (!state.currentCandidateSdf) {
    status.textContent = item.live_candidate_status ? tr("unavailable") : "";
    modelStatus.textContent = tr("emptyResult");
    return;
  }
  if (item.live_matches_candidate_snapshot === false) {
    status.textContent = tr("generatedSnapshotMismatch");
  } else if (item.live_matches_candidate_snapshot === true) {
    status.textContent = tr("generatedSnapshotMatch");
  } else {
    status.textContent = tr("generatedSnapshotIncomparable");
  }
  modelStatus.textContent = tr("rendered");
}

async function loadXyz(item, token) {
  const container = $("viewer3d");
  const technical = $("xyzTechnicalError");
  $("xyzText").textContent = tr("statusLoading");
  technical.hidden = true;
  technical.textContent = "";
  try {
    const xyz = await fetch(`/api/cases/${encodeURIComponent(item.case_id)}/xyz`).then((r) => {
      if (!r.ok) throw new Error(`XYZ load failed: ${r.status}`);
      return r.text();
    });
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    $("xyzText").textContent = xyz;
    state.xyzLoadStatus = "ok";
    state.xyzLoadError = "";
    render3d(xyz, token, item.case_id);
    renderReviewerSummary();
  } catch (error) {
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    state.xyzLoadStatus = "error";
    state.xyzLoadError = error.message;
    $("xyzText").textContent = error.message;
    container.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(tr("xyzUnavailable"))}</strong><span>${escapeHtml(tr("candidateUnavailableBecauseXyz"))}</span></div>`;
    technical.hidden = false;
    technical.textContent = msg("technicalDetail", { detail: error.message });
    renderReviewerSummary();
  }
}

async function loadCandidateSdf(item, token) {
  const text = $("candidateSdfText");
  const status = $("candidateSdfStatus");
  const modelStatus = $("candidateModelStatus");
  const technical = $("candidateTechnicalError");
  text.textContent = tr("statusLoading");
  status.textContent = "";
  modelStatus.textContent = "";
  state.currentCandidateSdf = "";
  state.currentLiveCandidate = null;
  technical.hidden = true;
  technical.textContent = "";
  try {
    const data = await api(
      `/api/cases/${encodeURIComponent(item.case_id)}/candidate-sdf?mode=${encodeURIComponent(state.twoDMode)}`,
    );
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    state.currentLiveCandidate = data;
    item.live_candidate_status = data.live_candidate_status || "";
    item.live_candidate_smiles = data.live_candidate_smiles || "";
    item.live_candidate_smiles_exact_match = data.live_candidate_smiles_exact_match;
    item.live_matches_candidate_snapshot = data.live_matches_candidate_snapshot;
    item.live_candidate_equivalence_method = data.live_candidate_equivalence_method || "";
    item.live_candidate_equivalence_reason = data.live_candidate_equivalence_reason || "";
    renderCaseHeader();
    renderReviewerSummary();
    renderDiagnostics();
    renderReviewerDetails();
    if (!data.available) {
      text.textContent = data.error || tr("reconstructionUnavailable");
      item.live_candidate_status = "unavailable";
      renderReviewerSummary();
      renderReviewerDetails();
      status.textContent = tr("unavailable");
      modelStatus.textContent = tr("unavailable");
      technical.hidden = !data.error;
      technical.textContent = data.error ? msg("technicalDetail", { detail: data.error }) : "";
      renderCandidate3d("", token, item.case_id);
      return;
    }
    state.currentCandidateSdf = data.sdf || "";
    text.textContent = state.currentCandidateSdf;
    if (data.live_matches_candidate_snapshot === false) {
      status.textContent = tr("generatedSnapshotMismatch");
    } else if (data.live_matches_candidate_snapshot === true) {
      status.textContent = tr("generatedSnapshotMatch");
    } else {
      status.textContent = state.currentCandidateSdf
        ? tr("generatedSnapshotIncomparable")
        : tr("emptyResult");
    }
    renderCandidate3d(state.currentCandidateSdf, token, item.case_id);
    modelStatus.textContent = state.currentCandidateSdf ? tr("rendered") : tr("emptyResult");
  } catch (error) {
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    item.live_candidate_status = "error";
    renderReviewerSummary();
    renderReviewerDetails();
    text.textContent = error.message;
    status.textContent = tr("error");
    modelStatus.textContent = tr("error");
    technical.hidden = false;
    technical.textContent = msg("technicalDetail", { detail: error.message });
    renderCandidate3d("", token, item.case_id);
  }
}

function mappingAmbiguityEvidence(data) {
  const locations = data.mapping_ambiguity_locations || {};
  const affected = Array.isArray(locations.affected_xyz_atoms)
    ? locations.affected_xyz_atoms.filter((index) => Number.isInteger(Number(index)))
    : [];
  const alternatives = Array.isArray(locations.alternatives) ? locations.alternatives : [];
  if (!locations.location_proven || !affected.length) {
    return {
      location: tr("ambiguityLocationUnknown"),
      alternatives: "",
    };
  }
  const formatDifference = (difference) => {
    if (difference.kind === "hydrogen_assignment") {
      const hydrogens = (difference.hydrogen_xyz_atoms || []).map((index) => `H/XYZ #${index}`).join(", ") || "H";
      const candidate = difference.candidate_center_xyz == null ? "—" : `XYZ #${difference.candidate_center_xyz}`;
      const reference = (difference.reference_center_xyz || []).map((index) => `XYZ #${index}`).join(", ") || "—";
      return `${hydrogens}: Candidate ${candidate} / Reference ${reference}`;
    }
    const atoms = (difference.xyz_atoms || []).map((index) => `XYZ #${index}`).join("–");
    const kind = difference.kind === "metal_bond" ? "metal coordination" : "organic bond";
    return `${kind} ${atoms}: Candidate ${difference.candidate_bond || "none"} / Reference ${difference.reference_bond || "none"}`;
  };
  const alternativeLines = alternatives
    .map((alternative) => {
      const descriptions = (alternative.differences || []).map(formatDifference).join("; ");
      return descriptions ? `${msg("mappingAlternative", { index: alternative.alternative })}: ${descriptions}` : "";
    })
    .filter(Boolean);
  return {
    location: affected.map((index) => `XYZ #${index}`).join(", "),
    alternatives: alternativeLines.join("\n"),
  };
}

async function loadReferenceXyz(item, token) {
  const container = $("referenceViewer3d");
  const technical = $("referenceXyzTechnicalError");
  state.referenceXyzSdf = "";
  state.referenceXyzFailure = "";
  technical.hidden = true;
  technical.textContent = "";
  technical.classList.remove("representative-mapping-warning");
  container.innerHTML = `<div class="empty">${escapeHtml(tr("statusLoading"))}</div>`;
  try {
    const data = await api(`/api/cases/${encodeURIComponent(item.case_id)}/reference-xyz`);
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    if (!data.available || !data.sdf) {
      const failure = String(data.failure_code || data.error || "atom_correspondence_not_reliable");
      state.referenceXyzFailure = data.error && data.error !== failure
        ? `${failure}: ${data.error}`
        : failure;
      const reviewerMessage = failure === "reference_missing_or_invalid"
        ? tr("referenceXyzMissing")
        : failure.includes("ambiguous") || failure === "atom_correspondence_not_reliable"
          ? tr("referenceXyzAmbiguous")
          : tr("referenceXyzUnreliable");
      const ambiguityEvidence = mappingAmbiguityEvidence(data);
      const ambiguity = data.mapping_confidence === "ambiguous"
        ? `<dl class="reference-xyz-ambiguity">
            <dt>${escapeHtml(tr("ambiguityType"))}</dt>
            <dd>${escapeHtml(tr(`mappingAmbiguityTypes.${data.mapping_ambiguity_type}`, data.mapping_ambiguity_type || "ambiguous"))}</dd>
            <dt>${escapeHtml(tr("unavailableReason"))}</dt>
            <dd>${escapeHtml(tr(`mappingAmbiguityReasons.${data.mapping_ambiguity_reason}`, reviewerMessage))}</dd>
            <dt>${escapeHtml(tr("ambiguityLocation"))}</dt>
            <dd><code>${escapeHtml(ambiguityEvidence.location)}</code></dd>
            ${ambiguityEvidence.alternatives ? `<dt>${escapeHtml(tr("ambiguityAlternatives"))}</dt><dd class="mapping-alternatives">${escapeHtml(ambiguityEvidence.alternatives)}</dd>` : ""}
          </dl>`
        : `<span>${escapeHtml(reviewerMessage)}</span>`;
      container.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(tr("referenceXyzUnavailable"))}</strong>${ambiguity}</div>`;
      technical.hidden = true;
      technical.textContent = "";
      state.referenceXyzViewer = null;
      renderDiagnostics();
      return;
    }
    state.referenceXyzSdf = data.sdf;
    state.mappedComparison = data.mapped_comparison || {};
    if (data.mapping_is_representative) {
      const ambiguityType = tr(
        `mappingAmbiguityTypes.${data.mapping_ambiguity_type}`,
        data.mapping_ambiguity_type || "ambiguous",
      );
      const ambiguityReason = tr(
        `mappingAmbiguityReasons.${data.mapping_ambiguity_reason}`,
        data.mapping_ambiguity_reason || tr("referenceXyzAmbiguous"),
      );
      technical.classList.add("representative-mapping-warning");
      technical.hidden = false;
      technical.innerHTML = `<strong>${escapeHtml(tr("representativeMappingWarning"))}</strong><span>${escapeHtml(tr("ambiguityType"))}: ${escapeHtml(ambiguityType)} · ${escapeHtml(tr("diagnosticReasonLabel"))}: ${escapeHtml(ambiguityReason)}</span>`;
    }
    renderMappedComparisonNote();
    renderReferenceXyz3d(data.sdf, token, item.case_id);
  } catch (error) {
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    container.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(tr("referenceXyzUnavailable"))}</strong><span>${escapeHtml(tr("referenceXyzUnreliable"))}</span></div>`;
    state.referenceXyzFailure = error.message;
    technical.hidden = true;
    technical.textContent = "";
    state.referenceXyzViewer = null;
    renderDiagnostics();
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

function mappedCoordinationEdges() {
  const edges = state.mappedComparison?.coordination_edges;
  return Array.isArray(edges)
    ? edges.filter((edge) => edge.presence !== "common" || (edge.mapped_ligand_group || []).length > 1)
    : [];
}

function renderMappedComparisonNote() {
  const note = $("mappedComparisonNote");
  if (state.xyzComparisonMode !== "mapped" || !mappedCoordinationEdges().length) {
    note.hidden = true;
    note.textContent = "";
    return;
  }
  const lines = mappedCoordinationEdges().map((edge) => {
    const pair = `${edge.metal_element} C/XYZ #${edge.candidate_metal_xyz_index} ↔ ${edge.donor_element} C/XYZ #${edge.candidate_donor_xyz_index}`;
    const distance = Number.isFinite(Number(edge.distance)) ? ` · ${Number(edge.distance).toFixed(3)} Å` : "";
    if (edge.presence === "common") {
      const group = edge.mapped_ligand_group || [];
      const correspondence = group.length > 1
        ? ` · ligand mapping ${group.map((atom) => `${atom.element} C/XYZ #${atom.candidate_xyz_index} ↔ R #${atom.reference_atom_index}`).join(", ")}`
        : "";
      return `${tr("mappedDonorPreserved")} · ${pair}${distance}${correspondence}`;
    }
    return `${edge.presence === "candidate_only" ? "Candidate-only" : "Reference-only"} coordination after mapped comparison · ${pair}${distance}`;
  });
  note.textContent = lines.join("\n");
  note.hidden = false;
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

function render3d(xyz, token = state.caseRequestToken, caseId = state.current?.case_id) {
  const container = $("viewer3d");
  container.innerHTML = "";
  if (!window.$3Dmol) {
    container.innerHTML = `<div class="empty">${escapeHtml(tr("threeDmolUnavailableXyz"))}</div>`;
    return;
  }
  requestAnimationFrame(() => {
    if (!isCurrentCaseRequest(token, caseId)) return;
    const viewer = window.$3Dmol.createViewer(container, { backgroundColor: "white" });
    const useMappedCandidate = state.xyzComparisonMode === "mapped" && state.currentCandidateSdf;
    viewer.addModel(useMappedCandidate ? state.currentCandidateSdf : xyz, useMappedCandidate ? "sdf" : "xyz");
    applyViewerStyle(viewer);
    viewer.zoomTo();
    if (typeof viewer.resize === "function") {
      viewer.resize();
    }
    viewer.render();
    state.xyzViewer = viewer;
    if (state.referenceXyzViewer) {
      copyViewerPose(viewer, state.referenceXyzViewer);
      state.referenceXyzViewer.render();
    }
    if (state.currentCandidateSdf) {
      renderCandidate3d(state.currentCandidateSdf, token, caseId);
    }
  });
}

function renderCandidate3d(sdf, token = state.caseRequestToken, caseId = state.current?.case_id) {
  const container = $("viewerCandidate3d");
  container.innerHTML = "";
  if (!sdf) {
    container.innerHTML = `<div class="empty">${escapeHtml(tr("emptyCandidate3d"))}</div>`;
    state.candidateViewer = null;
    return;
  }
  if (!window.$3Dmol) {
    container.innerHTML = `<div class="empty">${escapeHtml(tr("threeDmolUnavailableSdf"))}</div>`;
    state.candidateViewer = null;
    return;
  }
  requestAnimationFrame(() => {
    if (!isCurrentCaseRequest(token, caseId)) return;
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

function renderReferenceXyz3d(
  sdf,
  token = state.caseRequestToken,
  caseId = state.current?.case_id,
) {
  const container = $("referenceViewer3d");
  container.innerHTML = "";
  if (!sdf || !window.$3Dmol) {
    container.innerHTML = `<div class="empty">${escapeHtml(tr("referenceXyzUnavailable"))}</div>`;
    state.referenceXyzViewer = null;
    return;
  }
  requestAnimationFrame(() => {
    if (!isCurrentCaseRequest(token, caseId)) return;
    const viewer = window.$3Dmol.createViewer(container, { backgroundColor: "white" });
    viewer.addModel(sdf, "sdf");
    applyViewerStyle(viewer);
    viewer.zoomTo();
    copyViewerPose(state.xyzViewer, viewer);
    if (typeof viewer.resize === "function") viewer.resize();
    viewer.render();
    state.referenceXyzViewer = viewer;
  });
}

async function loadPair(token = state.caseRequestToken) {
  await Promise.all([
    loadRender("candidate", "primary", token),
    loadRender("reference", "secondary", token),
  ]);
}

function renderErrorDetail(kind, { httpStatus = null, payloadError = "", message = "" } = {}) {
  const parts = [`kind=${kind}`];
  if (httpStatus !== null && httpStatus !== undefined) parts.push(`HTTP ${httpStatus}`);
  if (payloadError) parts.push(`payload.error=${payloadError}`);
  if (message && message !== payloadError) parts.push(message);
  return parts.join(" · ");
}

function showReferenceMessage(message, technicalDetail = "") {
  const box = $("secondarySvg");
  box.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(message)}</strong></div>`;
  setImageZoomState(box, tr("reference"));
  const technical = $("referenceTechnicalError");
  technical.hidden = !technicalDetail;
  technical.textContent = technicalDetail
    ? msg("technicalDetail", { detail: technicalDetail })
    : "";
}

async function loadRender(kind, slot, token = state.caseRequestToken) {
  const item = state.current;
  if (!item || !isCurrentCaseRequest(token, item.case_id)) return;
  const title = slot === "primary" ? $("primaryRenderTitle") : $("secondaryRenderTitle");
  const reason =
    slot === "primary" ? $("primaryReferenceReason") : $("secondaryReferenceReason");
  const box = slot === "primary" ? $("primarySvg") : $("secondarySvg");
  const smiles = slot === "primary" ? $("primarySmiles") : $("secondarySmiles");
  title.textContent = tr(labels[kind], kind);
  const referenceReason = String(item?.reference_answer_reason || "").trim();
  reason.textContent = kind === "reference" && item?.reference_answer_wrong === "True"
    ? referenceReason
    : "";
  reason.hidden = !reason.textContent;
  box.innerHTML = `<div class="empty">${escapeHtml(tr("rendering"))}</div>`;
  setImageZoomState(box, title.textContent);
  smiles.textContent = "";
  if (kind === "reference") {
    $("referenceTechnicalError").hidden = true;
    $("referenceTechnicalError").textContent = "";
  }
  const initialReferenceState = kind === "reference" ? referenceState(item) : null;
  if (initialReferenceState === "missing" || initialReferenceState === "parse_invalid") {
    state.referenceRenderError = "";
    showReferenceMessage(
      initialReferenceState === "missing"
        ? tr("referenceMissingNotice")
        : tr("referenceInvalidNotice"),
    );
    renderReviewerSummary();
    return;
  }
  try {
    const data =
      kind === "candidate"
        ? state.currentLiveCandidate
        : await api(
            `/api/cases/${encodeURIComponent(item.case_id)}/render?kind=${encodeURIComponent(kind)}&mode=${encodeURIComponent(state.twoDMode)}&localize=1`,
          );
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    const renderError = kind === "candidate" ? data?.render_error || data?.error : data?.error;
    if (!data || renderError) {
      const detail = renderErrorDetail(kind, {
        payloadError: renderError || "",
        message: !data ? "empty response" : "",
      });
      if (kind === "candidate") {
        box.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(tr("candidateUnavailable"))}</strong></div>`;
        $("candidateTechnicalError").hidden = false;
        $("candidateTechnicalError").textContent = msg("technicalDetail", { detail });
      } else {
        state.referenceRenderStatus = "render_failed";
        state.referenceRenderError = detail;
        showReferenceMessage(tr("referenceRenderFailed"), detail);
      }
    } else {
      if (data.svg) {
        box.innerHTML = data.svg;
        if (kind === "reference") {
          state.referenceRenderStatus = "available";
          state.referenceRenderError = "";
          $("referenceTechnicalError").hidden = true;
          $("referenceTechnicalError").textContent = "";
        }
      } else {
        const detail = renderErrorDetail(kind, { message: "empty SVG" });
        if (kind === "reference") {
          state.referenceRenderStatus = "render_failed";
          state.referenceRenderError = detail;
          showReferenceMessage(tr("referenceRenderFailed"), detail);
        } else {
          box.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(tr("candidateUnavailable"))}</strong></div>`;
          $("candidateTechnicalError").hidden = false;
          $("candidateTechnicalError").textContent = msg("technicalDetail", { detail });
        }
      }
    }
    setImageZoomState(box, title.textContent);
    smiles.textContent = kind === "candidate" ? data?.live_candidate_smiles || "" : data?.smiles || "";
    if (kind === "reference") {
      renderCaseHeader();
      renderReviewerSummary();
    }
  } catch (error) {
    if (!isCurrentCaseRequest(token, item.case_id)) return;
    const detail = renderErrorDetail(kind, {
      httpStatus: error.httpStatus,
      payloadError: error.payloadError,
      message: error.message,
    });
    box.innerHTML = `<div class="reviewer-error"><strong>${escapeHtml(kind === "candidate" ? tr("candidateUnavailable") : tr("referenceRenderFailed"))}</strong></div>`;
    setImageZoomState(box, title.textContent);
    if (kind === "reference") {
      state.referenceRenderStatus = "render_failed";
      state.referenceRenderError = detail;
      showReferenceMessage(tr("referenceRenderFailed"), detail);
      renderCaseHeader();
      renderReviewerSummary();
    } else {
      $("candidateTechnicalError").hidden = false;
      $("candidateTechnicalError").textContent = msg("technicalDetail", { detail });
    }
  }
}

function renderReviewHistory() {
  const list = $("recentReviewList");
  const active = state.reviewHistory.find((item) => !item.undone);
  $("undoReview").disabled = !active || state.savingReview || state.undoingReview;
  $("lastReviewSummary").textContent = active
    ? (active.mutation_type === "family_qa"
      ? `${active.family_id} → ${active.status}`
      : msg("justReviewed", { caseId: active.case_id, status: active.status }))
    : "";
  if (!state.reviewHistory.length) {
    list.innerHTML = `<p class="muted">${escapeHtml(tr("noRecentReviews"))}</p>`;
    return;
  }
  list.innerHTML = state.reviewHistory.map((item) => `
    <button class="recent-review-item${item.undone ? " undone" : ""}" data-case-id="${escapeHtml(item.case_id)}" type="button">
      <strong>${escapeHtml(item.case_id)}</strong>
      <span>${escapeHtml(item.status)}</span>
      <span>${escapeHtml(item.reviewer || "—")}</span>
      <time>${escapeHtml(new Date(item.timestamp).toLocaleTimeString())}</time>
    </button>
  `).join("");
  list.querySelectorAll(".recent-review-item").forEach((button) => {
    button.addEventListener("click", () => loadCase(button.dataset.caseId));
  });
}

async function loadReviewHistory() {
  const data = await api("/api/review-history");
  state.reviewHistory = data.items || [];
  renderReviewHistory();
}

async function undoLastReview() {
  if (state.undoingReview || state.savingReview) return;
  const latest = state.reviewHistory.find((item) => !item.undone);
  if (!latest) return;
  state.undoingReview = true;
  renderReviewHistory();
  try {
    const result = await api("/api/review-undo", {
      method: "POST",
      body: JSON.stringify({ mutation_id: latest.mutation_id }),
    });
    $("saveState").textContent = msg("undoComplete", { caseId: result.case_id });
    await Promise.all([loadStats(), loadReviewReasons(), loadReviewHistory(), loadFamilyQa()]);
    await loadCases();
    if (result.case_id) await loadCase(result.case_id);
    else {
      const targetCase = activeFamily()?.representatives[0]?.case_id;
      if (targetCase) await loadCase(targetCase);
    }
  } catch (error) {
    $("saveState").textContent = error.message;
    await loadReviewHistory();
  } finally {
    state.undoingReview = false;
    renderReviewHistory();
  }
}

async function saveReview(status) {
  if (!state.current || state.savingReview || activeFamily()) return;
  state.savingReview = true;
  const currentIndex = state.cases.findIndex((item) => item.case_id === state.current.case_id);
  const queuedNextCaseId = currentIndex >= 0 ? state.cases[currentIndex + 1]?.case_id : "";
  $("saveState").textContent = tr("saving");
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
    const result = await api(`/api/cases/${encodeURIComponent(state.current.case_id)}/review`, {
      method: "POST",
      body: JSON.stringify(payload),
    });
    if (result.mutation) {
      state.reviewHistory.unshift(result.mutation);
      state.reviewHistory = state.reviewHistory.slice(0, 20);
      renderReviewHistory();
    }
    $("saveState").textContent = result.fixture
      ? msg("savedFixture", { file: result.fixture.structure_file })
      : tr("savedNoFixture");
    await Promise.all([loadStats(), loadReviewReasons()]);
    await loadCases();
    if (queuedNextCaseId) {
      await loadCase(queuedNextCaseId);
      return;
    }
    const nextIndex = currentIndex >= 0 ? Math.min(currentIndex, state.cases.length - 1) : 0;
    const nextCase = state.cases[nextIndex];
    if (nextCase) await loadCase(nextCase.case_id);
  } catch (error) {
    $("saveState").textContent = error.message;
  } finally {
    state.savingReview = false;
  }
}

async function navigateCase(delta) {
  if (!state.cases.length) return;
  const currentIndex = state.current
    ? state.cases.findIndex((item) => item.case_id === state.current.case_id)
    : -1;
  const targetIndex = currentIndex + delta;
  if (targetIndex >= 0 && targetIndex < state.cases.length) {
    await loadCase(state.cases[targetIndex].case_id);
    return;
  }
  if (delta > 0 && state.offset + state.limit < state.total) {
    state.offset += state.limit;
    await loadCases();
    if (state.cases.length) await loadCase(state.cases[0].case_id);
  } else if (delta < 0 && state.offset > 0) {
    state.offset = Math.max(0, state.offset - state.limit);
    await loadCases();
    if (state.cases.length) await loadCase(state.cases[state.cases.length - 1].case_id);
  }
}

function reviewShortcut(event) {
  const target = event.target;
  if (target?.matches?.("input, textarea, select, [contenteditable='true']")) return;
  if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === "z") {
    event.preventDefault();
    undoLastReview();
    return;
  }
  const statusByKey = {
    1: "accept_candidate",
    2: "accept_reference",
    3: "accept_both",
    4: "manual_reference",
    5: "reference_answer_wrong",
    6: "needs_followup",
    7: "skip",
  };
  if (statusByKey[event.key]) {
    if (activeFamily()) return;
    event.preventDefault();
    saveReview(statusByKey[event.key]);
    return;
  }
  if (event.key === "ArrowLeft") {
    event.preventDefault();
    navigateCase(-1);
  } else if (event.key === "ArrowRight") {
    event.preventDefault();
    navigateCase(1);
  }
}

async function removeCurrentFixture() {
  const item = state.current;
  if (!item?.fixture) return;
  const confirmed = window.confirm(
    msg("removeFixtureConfirm", { file: item.fixture.structure_file, caseId: item.case_id }),
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
      setKetcherStatus("ketcherReady");
      return ketcher;
    }
    await new Promise((resolve) => setTimeout(resolve, 250));
  }
  throw localizedError("ketcherNotReady");
}

async function setKetcherMolecule(smiles) {
  if (!smiles) {
    setKetcherStatus("noLoadableSmiles");
    return;
  }
  setKetcherStatus("loadingMolecule");
  try {
    const ketcher = await waitForKetcher();
    if (typeof ketcher.setMolecule !== "function") {
      throw localizedError("ketcherNoSetMolecule");
    }
    await ketcher.setMolecule(smiles);
    setKetcherStatus("moleculeLoaded");
  } catch (error) {
    setKetcherStatus(error.i18nKey || "", {}, error.message);
  }
}

async function readKetcherToForm({ rethrow = false, requireSmiles = false } = {}) {
  setKetcherStatus("readingCanvas");
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
    if (requireSmiles && !smiles) throw localizedError("manualSmilesRequired");
    if (!molblock.trim() && !smiles) {
      if (molblockError) throw molblockError;
      throw localizedError("canvasEmpty");
    }
    setKetcherStatus("canvasCopied");
    return { molblock, smiles };
  } catch (error) {
    setKetcherStatus(error.i18nKey || "", {}, error.message);
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
  applyLanguage();
  window.addEventListener("load", () => {
    if (window.$3Dmol && $("xyzText").textContent) {
      render3d($("xyzText").textContent);
    }
    if (window.$3Dmol && state.referenceXyzSdf) {
      renderReferenceXyz3d(state.referenceXyzSdf);
    }
  });
  $("ketcherFrame").addEventListener("load", () => {
    waitForKetcher(30000)
      .then(() => {
        const workspace = document.querySelector(".workspace");
        if (workspace) workspace.scrollTop = 0;
      })
      .catch((error) => {
        setKetcherStatus(error.i18nKey || "", {}, error.message);
      });
  });
  $("refreshStats").addEventListener("click", () => {
    loadStats();
    loadReviewReasons();
    loadCases(true);
  });
  $("languageToggle").addEventListener("click", toggleLanguage);
  $("reviewer").addEventListener("focus", () => loadReviewReasons());
  $("reviewerDetails").addEventListener("toggle", () => {
    if ($("reviewerDetails").open) loadGraphEvidence(state.caseRequestToken);
  });
  $("jumpToReviewerDetails").addEventListener("click", jumpToReviewerDetails);
  document.querySelectorAll(".show-full-graph").forEach((checkbox) => {
    checkbox.addEventListener("change", () => renderGraphInspector(checkbox.dataset.kind));
  });
  $("categoryFilter").addEventListener("change", () => loadCases(true));
  $("statusFilter").addEventListener("change", () => loadCases(true));
  $("reviewReasonFilter").addEventListener("change", () => {
    if ($("reviewReasonFilter").value && $("statusFilter").value === "unreviewed") {
      $("statusFilter").value = "";
    }
    loadCases(true);
  });
  $("triageFilter").addEventListener("change", () => loadCases(true));
  $("familyQueue").addEventListener("change", async () => {
    state.activeFamilyId = $("familyQueue").value;
    await loadCases(true);
    if (state.cases.length) await loadCase(state.cases[0].case_id);
  });
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
  $("undoReview").addEventListener("click", undoLastReview);
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
      if (state.referenceXyzSdf) renderReferenceXyz3d(state.referenceXyzSdf);
      if (state.currentCandidateSdf) renderCandidate3d(state.currentCandidateSdf);
    });
  });
  document.querySelectorAll(".xyz-comparison-mode").forEach((button) => {
    button.addEventListener("click", () => {
      if (button.dataset.mode === state.xyzComparisonMode) return;
      state.xyzComparisonMode = button.dataset.mode;
      document.querySelectorAll(".xyz-comparison-mode").forEach((item) => {
        item.classList.toggle("active", item.dataset.mode === state.xyzComparisonMode);
      });
      renderMappedComparisonNote();
      if ($("xyzText").textContent) render3d($("xyzText").textContent);
      if (state.referenceXyzSdf) renderReferenceXyz3d(state.referenceXyzSdf);
    });
  });
  document.querySelectorAll(".two-d-mode").forEach((button) => {
    button.addEventListener("click", async () => {
      if (!state.current || button.dataset.mode === state.twoDMode) return;
      document.querySelectorAll(".two-d-mode").forEach((item) => item.classList.remove("active"));
      button.classList.add("active");
      state.twoDMode = button.dataset.mode;
      const token = state.caseRequestToken;
      await loadCandidateSdf(state.current, token);
      if (isCurrentCaseRequest(token, state.current.case_id)) await loadPair(token);
    });
  });
  document.querySelectorAll(".decision").forEach((button) => {
    button.addEventListener("click", () => saveReview(button.dataset.status));
  });
  document.querySelectorAll("[data-rep-mark]").forEach((button) => {
    button.addEventListener("click", () => mutateFamilyQa("representative_mark", button.dataset.repMark));
  });
  document.querySelectorAll("[data-family-decision]").forEach((button) => {
    button.addEventListener("click", () => mutateFamilyQa("decision", button.dataset.familyDecision));
  });
  document.addEventListener("keydown", reviewShortcut);
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
  await loadReviewReasons();
  await loadReviewHistory();
  await loadFamilyQa({ preserveFamily: false });
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
