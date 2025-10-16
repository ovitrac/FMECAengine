# 🛡️ FMECAengine — Failure Mode, Effects & Criticality Analysis for Mass Transfer

[![DeepWiki – FMECAengine](https://img.shields.io/badge/DeepWiki-FMECAengine-0969da?style=flat&logo=bookstack&logoColor=white)](https://deepwiki.com/ovitrac/FMECAengine)
[![License: CeCILL-B](https://img.shields.io/badge/License-CeCILL--B-blue.svg)](#-license)
[![MATLAB/Octave](https://img.shields.io/badge/MATLAB%2FOctave-supported-success.svg)](#-requirements)

**FMECAengine** is a scientific software toolkit to perform **FMECA** in **mass-transfer/migration** problems (e.g., food-contact materials, porous media, barrier design, supply-chain scenarios).  
It unifies *failure-mode identification*, *severity/occurrence/detectability scoring*, and *physics-based transport modeling*—with an **internal inference language** `key2key()` enabling **rule-enforced AI workflows** and **rapid orchestration of complex scenarios**.

---

## ✨ What’s inside

- **Domain-adapted FMECA** for transfer/migration pathways (materials, steps, environments).
- **Criticality ranking** (S×O×D or custom functions), thresholds, and priority lists.
- **Coupling to transport kernels** (diffusion/partition models, kinetics).
- **Scenario automation via `key2key()`**: a compact language to query/join heterogeneous tables, propagate constraints, evaluate formulas, and emit a *query graph* for tracing/explainability.
- **Batch & Monte-Carlo style exploration**: run thousands of configurations or entire supply-chain variants using concise keys.
- **Tabular/JSON exports** and figures for reports.

---

## 🧠 The `key2key()` language — rule-enforced orchestration for GenAI

`key2key()` is a **multivalued, relational inference language** embedded in FMECAengine (see `key2key.m`).  
It composes *table-to-table relationships* (one-to-one, one-to-many, many-to-many, cascades) with an intuitive syntax:

```
<selector>:<tableA>::<colA>-><tableB>::<colB>-><tableC>::<colC> ...
```

### Key properties that make it valuable for AI agents

- **Deterministic joins with built-in uniqueness** on string-valued chains to avoid ambiguous expansions.
- **Lists, alternations, and regex** directly in keys (`|`, `\d`, character classes, quantifiers with `{{}}`, etc.).
- **Mathematical expressions** embedded in keys (e.g., `min(...)`, `max(...)`, calling domain models such as `Dpiringer(...)`).
- **Chained/cascaded joins** across multiple tables and *double/triple keys* for scenario coupling (e.g., conditions, materials, and contact states).
- **Explainability**: returns a **query tree** (`key2keygraph`) for traceable, human-auditable pipelines (ideal for GenAI *tool-use* and *chain-of-thought externalization*).

### Minimal examples (see `key2key.m` for many more)

```matlab
% Substances used with a polymer class, then fetch molar masses M:
'PP:polymer::name->classadditives:substance::class->name:substance::name->M'

% Alternatives & regex:
'\antiUV\d\ :substance::name->M'         % any antiUV<digit>
'antiUV2 | antiox1 :substance::name->M'  % union

% Invoke a physical model inside the key:
'Dpiringer(LDPE, \antiUV{3|2}\ :substance::name->M, 40)'

% Double-key example: polymer comes from a scenario table, T from a contact table
'Dpiringer((scenarioA:scenario::id->polymer), (PP:polymer::name->...->M), (cond3:contact::condition->temperature))'
```

> 🔎 **Why this matters for GenAI**
> `key2key()` gives LLM agents a *formal, minimal* instruction set to **enforce constraints, propagate rules**, and **compose multi-table joins and physics calls** safely—without free-form SQL. It is thus well-suited for **agentic pipelines** that must remain **deterministic, auditable, and reproducible**.

------

## 📦 Install / Run

> Works with **MATLAB** and **GNU Octave**.

```matlab
git clone https://github.com/ovitrac/FMECAengine.git
cd FMECAengine
# (MATLAB/Octave) addpath(pwd), then run examples or your own scripts
```


Typical entry points:

- `fmecaengine.m`, `fmecasingle.m`, `fmecagraph.m`, `fmecamerge.m`
- `key2key.m` + `key2keygraph.m` (inference & explainability)
- Diffusion/partition helpers: `Dpiringer.m`, `Dhelmroth.m`, `Dfuller.m`, …

------

## 🧩 Quick workflow

1. **Describe** components/materials/steps and potential *failure (transfer) modes*.
2. **Score** each mode: severity **S**, occurrence **O**, detectability **D** (or custom).
3. **Quantify** transport via embedded models (diffusion, partition, kinetics).
4. **Rank** by criticality (e.g., $C = S \times O \times D$) and apply thresholds.
5. **Automate scenarios** with `key2key()` across tables (materials ↔ substances ↔ conditions), including regex and lists.
6. **Export and review** the ranked modes and traces (tables/JSON/graphs).

------

## 📚 Publications using **FMECAengine** (Vitrac & Nguyen)

> *This list is curated from public records (Google Scholar/Scopus-indexed pages) and references explicitly citing or describing **FMECAengine** / **key2key**. Please open a PR to extend.*

- **Nguyen, P-M., Goujon, A., Sauvegrain, P., Vitrac, O.** (2013). *A computer-aided methodology to design safe food packaging and related systems.* **AIChE Journal**, 59(4), 1288–1306. doi: **10.1002/aic.14056**.
  (Cites the **FMECAengine** open-source project.)
- **Zhu, Y., Vitrac, O., Nguyen, P-M., Hayert, M.** (2019). *Toward Safer and Ecodesigned Food Packaging Systems.* **Frontiers in Chemistry**, 7:349. doi: **10.3389/fchem.2019.00349**.
  (Explicit GitHub citation of **FMECAengine**.)
- **Nguyen, P-M., Zhu, Y., Vitrac, O.** (2019). *The Ubiquitous Issue of Cross-Mass Transfer: Applications to Single-Use Systems.* (Open-access article referencing **FMECAengine**.)
- **Vitrac, O.** (2014). *Food Packaging: New Directions for the Control of Additive Migration.* In **Wiley** book chapter.
  (Includes a direct reference to **FMECAengine**.)
- **Nguyen, P-M., Berrard, C., Daoud, N., Vitrac, O.** (2024). *Assessment of chemical risks and circular economy implications of recycled PET in food packaging with functional barriers.* **Cleaner Materials** (Elsevier).
  (Mentions **FMECAengine 0.63** as FMECA applied to mass transfer.)
- **Vitrac, O., Nguyen, P-M., Hayert, M.** (2022). *In Silico Prediction of Food Properties: A Multiscale Perspective.* **Frontiers in Chemical Engineering**, 3:786879.
  (Broader multiscale framework citing **FMECAengine / key2key** concepts and 3D extensions.)

> Additional institutional and training materials by **INRAE** and **FitNESS** reference FMECAengine and `key2key()` for safe-by-design packaging and training programs.

------

## 🗺️ Where `key2key()` shines (typical patterns)

- **Functional-barrier design**: connect *polymer → additive classes → substances → properties (M, D, k)*; invoke `Dpiringer(...)`/others to compute permeation bounds under *T, t, geometry*.
- **Supply-chain composition**: material → process step → environment → distribution; run thousands of *what-ifs* by varying tables and using alternations/regex.
- **Explainable AI**: emit a *query tree* and visualize joins via `key2keygraph`, enabling human review and audit.

------

## 🧪 Minimal snippet

```matlab
addpath(genpath(pwd));                         % one-shot load
db = loadFMECAEngineDB('mydb.ods');            % or assemble struct tables programmatically

% Fetch M for all anti-UV substances used with PP and compute D at 40°C:
key = 'Dpiringer(PP, (PP:polymer::name->classadditives:substance::class->name:substance::name->M), 40)';
[val, ~, tree] = key2key(db, key);
key2keygraph(tree);                            % explain joins + values
```

------

## 🧷 Requirements

- **MATLAB** R2013b+ or **GNU Octave** (tested on multiple releases).
- Standard toolboxes for plotting/IO (see repo for helper scripts).

------

## 🤝 Contributing

1. Use issues/PRs for bugs and feature requests.
2. For new physical sub-models (diffusion/partition/kinetics), include **docstrings**, **units**, **domain of validity**, and **tests**.
3. For `key2key()` extensions, document **operators**, **regex coverage**, **uniqueness rules**, and **graph outputs**.

------

## 📄 License

**CeCILL-B** (GPL-compatible, French law). See `LICENSE`.

© 2011–2025 Olivier Vitrac and contributors.

