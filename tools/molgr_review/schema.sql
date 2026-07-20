PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS cases (
    case_id TEXT PRIMARY KEY,
    row_index INTEGER,
    source TEXT NOT NULL DEFAULT '',
    category TEXT NOT NULL DEFAULT '',
    xyz_path TEXT NOT NULL,
    total_charge INTEGER NOT NULL DEFAULT 0,
    total_radical_electrons INTEGER NOT NULL DEFAULT 0,
    spin_multiplicity INTEGER NOT NULL DEFAULT 1,
    reference_smiles TEXT NOT NULL DEFAULT '',
    candidate_smiles TEXT NOT NULL DEFAULT '',
    candidate_organic_smiles TEXT NOT NULL DEFAULT '',
    reference_organic_smiles TEXT NOT NULL DEFAULT '',
    candidate_status TEXT NOT NULL DEFAULT '',
    metadata_json TEXT NOT NULL DEFAULT '{}'
);

CREATE TABLE IF NOT EXISTS reviews (
    case_id TEXT PRIMARY KEY REFERENCES cases(case_id) ON DELETE CASCADE,
    status TEXT NOT NULL,
    corrected_smiles TEXT NOT NULL DEFAULT '',
    corrected_molblock TEXT NOT NULL DEFAULT '',
    notes TEXT NOT NULL DEFAULT '',
    reviewer TEXT NOT NULL DEFAULT '',
    updated_at TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS render_cache (
    case_id TEXT NOT NULL REFERENCES cases(case_id) ON DELETE CASCADE,
    kind TEXT NOT NULL,
    svg TEXT NOT NULL DEFAULT '',
    smiles TEXT NOT NULL DEFAULT '',
    error TEXT NOT NULL DEFAULT '',
    generated_at TEXT NOT NULL,
    PRIMARY KEY(case_id, kind)
);

CREATE INDEX IF NOT EXISTS idx_cases_category ON cases(category);
CREATE INDEX IF NOT EXISTS idx_reviews_status ON reviews(status);
