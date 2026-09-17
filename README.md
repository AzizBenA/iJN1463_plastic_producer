# iJN1463: metabolic modeling of mixed-plastic upcycling

Genome-scale metabolic models and flux balance analysis of plastic-derived monomer assimilation and (R)-3-hydroxybutyrate (R-3HB) production in *Pseudomonas putida* KT2440.

## Associated publication

This work was published as a **bioRxiv preprint**:

> Meng, H., et al. (2026). **Engineering Pseudomonas putida KT2440 for open-loop upcycling of mixed plastics.** *bioRxiv*, version 1, posted 25 March 2026. [doi:10.64898/2026.03.23.713816](https://doi.org/10.64898/2026.03.23.713816).

The [preprint](https://www.biorxiv.org/content/10.64898/2026.03.23.713816v1.abstract) identifies this repository as the source of the models and code used in its computational analysis. It describes engineering *P. putida* for plastic-monomer utilization and R-3HB production. The repository provides the metabolic modeling component; experimental methods and strain characterization are described in the manuscript.

## Scientific scope

The models extend the iJN1463 reconstruction with pathways connecting polyester-derived monomers to central carbon metabolism. The analysis notebook uses flux balance analysis (FBA) to maximize R-3HB production and examine how monomer composition affects theoretical yield.

| Abbreviation | Monomer | Uptake reaction in the notebook | Carbon atoms per molecule |
| --- | --- | --- | ---: |
| EG | Ethylene glycol | `SK_etheglycol_c` | 2 |
| TA | Terephthalic acid | `SK_terepa_c` | 8 |
| AA | Adipic acid | `SK_adpac_c` | 6 |
| BDO | 1,4-Butanediol | `SK_14btdl_c` | 4 |

Mixture simulations combine these monomers using compositions assigned to polyethylene terephthalate (PET), polyester-polyurethane (PU) soft segments, and poly(butylene adipate-co-terephthalate) (PBAT). These compositions are modeling assumptions for the represented materials, rather than universal polymer formulations. The notebook focuses on the four monomers listed above.

## Repository organization

| Resource | Contents |
| --- | --- |
| [Models/](Models/) | Six dated SBML model snapshots, from `updated_1` to `updated_6`. |
| [Scripts/Hao_yield_calculation.ipynb](Scripts/Hao_yield_calculation.ipynb) | Interactive model inspection, constraint changes, yield calculations, and polymer-mixture sweeps. |
| [Scripts/jupyter_utils.py](Scripts/jupyter_utils.py) | Helper functions for model inspection and tabular exports. |
| [Data/250326_yield_Data.csv](Data/250326_yield_Data.csv) | Historical single-monomer growth, yield, and uptake summary. |
| [Results/251218_polymer_mixture_yield_sweep_half.xlsx](Results/251218_polymer_mixture_yield_sweep_half.xlsx) | Archived polymer-mixture sweep results. |
| [Methods and reproducibility](docs/methods.md) | Notebook workflow, yield definitions, and interpretation of saved outputs. |
| [Model history](docs/model-history.md) | Snapshot provenance, reaction changes, and model-specific constraints. |
| [CITATION.bib](CITATION.bib) | Bibliographic record for the associated preprint. |

## Getting started

### Select a model

The notebook loads [250709_iJN1463_updated_5.sbml](Models/250709_iJN1463_updated_5.sbml). The newest archived snapshot is [260223_iJN1463_updated_6.sbml](Models/260223_iJN1463_updated_6.sbml). These snapshots have different reaction bounds; specify the exact filename when reporting an analysis. See the [model history](docs/model-history.md) before substituting one for the other.

### Prepare the environment

The notebook records Python **3.11.7**. The following packages cover its imports and Excel exports; the repository does not include a pinned environment or dependency lockfile.

```bash
python -m venv .venv
# Activate the environment before installing packages:
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install cobra pandas numpy requests openpyxl jupyterlab ipykernel
```

### Open the notebook

From the repository root, launch Jupyter in `Scripts/` so that the local helper import and relative model paths resolve correctly:

```bash
cd Scripts
python -m jupyterlab Hao_yield_calculation.ipynb
```

Use a kernel from the environment created above. Before executing the notebook, follow the [execution guide](docs/methods.md#notebook-execution): an intermediate cell writes directly to the archived `updated_6` model, and later cells represent separate historical sweep implementations.

To inspect a model independently, run this Python example from the repository root:

```python
from pathlib import Path
from cobra.io import read_sbml_model

model_path = Path("Models") / "260223_iJN1463_updated_6.sbml"
model = read_sbml_model(str(model_path))
print(model.id)
print(f"{len(model.reactions)} reactions; {len(model.metabolites)} metabolites")
print(model.objective.expression)
```

## Interpretation and reproducibility

FBA results describe feasible steady-state fluxes under the selected objective and reaction bounds. Product-maximizing solutions are theoretical predictions and should be interpreted separately from experimental titers, growth dynamics, and substrate-consumption order.

The notebook is an interactive research record. Its model state depends on execution order, and some historical output labels differ from the quantities calculated. The [methods guide](docs/methods.md) documents these details, including the distinction between carbon yield and molar product yield. For reuse, record the repository revision, SBML filename, solver and package versions, objective, uptake limits, and any additional constraints.

## Citation

Please cite the [associated preprint](https://doi.org/10.64898/2026.03.23.713816) when using these models or analyses. Import [CITATION.bib](CITATION.bib) into a reference manager for the full author list and publication metadata. Report the repository revision and model snapshot alongside the citation to identify the computational materials used.
