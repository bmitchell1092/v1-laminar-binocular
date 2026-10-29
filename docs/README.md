# Documentation

Laminar V1 recordings (macaques B and J, 2021–2022) during binocular,
monocular and dichoptic gratings at varying contrast. Also flash (evp), RF
mapping, orientation tuning and darkrest. Blackrock NSP + MonkeyLogic 2.

Started 2026-10-28 from the legacy [BM-ephys](https://github.com/bmitchell1092/BM-ephys)
documentation (tag `legacy-final`). Docs 01, 03 and 04 describe **legacy code**
and are kept for reference. Docs 02, 05 and 06 are maintained here.

| doc | contents | status |
|---|---|---|
| [01_preprocessing_pipeline.md](01_preprocessing_pipeline.md) | legacy raw → STIM → neural data: event codes, trial matching, photodiode, eye/contrast coding, V1 limits, NEV↔NSx clock | legacy reference |
| [02_dataset.md](02_dataset.md) | `D:\NEUROPHYS DATA\2021`: subjects, sessions, paradigms, channel maps, dual-probe days, processed folders | maintained |
| [03_project_history_and_conventions.md](03_project_history_and_conventions.md) | iScience 2022 → OD → LMN; legacy conventions, signal recipes, RF/orientation methods | legacy reference |
| [04_LMN_analysis.md](04_LMN_analysis.md) | legacy laminar-normalization analysis chain and figure list | legacy reference |
| [05_known_issues_and_QC.md](05_known_issues_and_QC.md) | issue register and verification plan | **maintained** |
| [06_rebuild_plan.md](06_rebuild_plan.md) | phases, decisions, to-dos | **maintained** |
| [07_supersaturation_assessment.md](07_supersaturation_assessment.md) | honest assessment: supersaturation vs depth/compartment (2026-10-28) |
| [08_binocular_feasibility.md](08_binocular_feasibility.md) | first-pass feasibility of the top literature questions (2026-10-29) |
| [../reports/V1 binocular contrast normalization.md](../reports/V1%20binocular%20contrast%20normalization.md) | literature sweep (2026-10-29): binocular/monocular contrast, interocular normalization, laminar timing; 16 ranked open questions vs this dataset (notes in `research_notes/`) |
| [legacy_TuneList_all_versions.csv](legacy_TuneList_all_versions.csv) | recovered legacy session metadata, every historical version | legacy reference |
