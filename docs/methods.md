# Methods and reproducibility

[Return to the repository overview](../README.md) · [Model history](model-history.md)

## Analysis materials

The analysis is implemented in [Hao_yield_calculation.ipynb](../Scripts/Hao_yield_calculation.ipynb), with supporting functions in [jupyter_utils.py](../Scripts/jupyter_utils.py). The notebook loads [updated_5](../Models/250709_iJN1463_updated_5.sbml) and modifies the model in memory. The [updated_6 snapshot](../Models/260223_iJN1463_updated_6.sbml) is also archived, but is not the notebook's input.

The [associated preprint](https://www.biorxiv.org/content/10.64898/2026.03.23.713816v1.abstract) provides the biological context and experimental methods. The details below describe the files and calculations currently present in this repository.

## Computational formulation

The notebook uses COBRApy to solve a steady-state flux balance problem, maximizing the flux through the product sink `SK_bhb_c`:

$$
\max_v\; v_{\mathrm{SK\_bhb\_c}}
\quad \text{subject to}\quad Sv = 0,\quad l \leq v \leq u.
$$

Here, $S$ is the stoichiometric matrix, $v$ is the reaction-flux vector, and $l$ and $u$ are the lower and upper bounds. The product metabolite is `bhb_c`, representing R-3HB. No positive minimum biomass flux is imposed in the notebook's product-yield calculations.

The four monomers enter through intracellular sink reactions. Negative sink flux represents uptake; the positive consumption rate is $u_i = -v_i$. Bounds of `(-capacity, 0)` specify an uptake ceiling, so the solver can consume less than the supplied capacity. They do not enforce exact uptake ratios.

### Initial substrate conditions

The initial notebook cell sets the following limits, in mmol gDW⁻¹ h⁻¹, where gDW denotes grams of cell dry weight:

| Reaction | Monomer | Lower bound | Upper bound |
| --- | --- | ---: | ---: |
| `SK_etheglycol_c` | EG | -3.3175 | 0 |
| `SK_terepa_c` | TA | -1.6500 | 0 |
| `SK_adpac_c` | AA | -3.3350 | 0 |
| `SK_14btdl_c` | BDO | -1.6675 | 0 |

That cell also closes `EX_h2s_e`. The input model has oxygen bounds `EX_o2_e = (-100, 0)`, a closed glucose exchange, and ATP maintenance `ATPM = (0.92, 0.92)`. Later mixture cells replace the four monomer uptake limits.

A separate, later constraint cell sets `SUCOAS = (-1000, 0)`, `PPC = (0, 1000)`, `PC = (0, 1000)`, and `THRA = (0, 0)`. In particular, this opens `PC` and closes `THRA` relative to the input snapshot. Early yield cells precede these changes. See the [snapshot bounds comparison](model-history.md) when selecting the conditions for an analysis.

## Polymer-mixture representation

The sweep functions combine monomer compositions as follows:

| Represented feedstock | EG | TA | AA | BDO |
| --- | ---: | ---: | ---: | ---: |
| PET | 0.500 | 0.500 | 0.000 | 0.000 |
| PU soft segments | 0.250 | 0.000 | 0.500 | 0.250 |
| PBAT | 0.000 | 0.222 | 0.278 | 0.500 |

Entries are monomer molar fractions within each modeled feedstock. Mixture weights allocate the total monomer uptake capacity across these composition vectors; they are not polymer mass fractions.

For monomer $i$, the capacity is $U_i = U_{\mathrm{total}} \sum_p f_p a_{pi}$, where $f_p$ is the mixture weight for polymer $p$ and $a_{pi}$ is its monomer fraction. The default total capacity is **10 mmol gDW⁻¹ h⁻¹**.

In the final sweep, each polymer is fixed in turn at a mixture weight of **0.5**. One remaining weight varies from **0 to 0.5**, inclusive, in increments of **0.01**, and the third makes the weights sum to one. This produces **153 scenarios**. An earlier sweep cell excludes the upper endpoint and produces **150 scenarios**.

The notebook also adds `PET_c`, `PU_c`, and `PBAT_c` pseudo-metabolites and corresponding assembly reactions. The mixture sweeps supply monomers through their sink bounds; the added assembly reactions are not an explicit model of enzymatic depolymerization.

## Yield definitions and units

Let $v_P$ denote the R-3HB sink flux, $u_i$ the actual positive consumption rate of monomer $i$, and $c_i$ its number of carbon atoms. R-3HB contains four carbon atoms.

| Quantity | Definition | Units |
| --- | --- | --- |
| Product flux | $v_P$ | mmol R-3HB gDW⁻¹ h⁻¹ |
| Carbon yield | $Y_C = 4v_P / \sum_i c_i u_i$ | mol product carbon per mol consumed substrate carbon |
| Molar product yield | $Y_{P/S} = v_P / \sum_i u_i$ | mol R-3HB per mol consumed monomer |

The notebook's `compute_carbon_yield` function uses the carbon-yield definition above. Its `compute_mol_yield` function calculates **$4v_P / \sum_i u_i$**, retaining the product-carbon multiplier. Consequently, the current function's exported `mol/mol_yield` column represents **mol product carbon per mol consumed monomer**. Dividing that calculated quantity by four gives the molar R-3HB yield. This distinction matters when comparing notebook exports with publication figures or experimental measurements.

The earlier cell that divides product carbon only by adipate carbon is specific to an adipate-only interpretation. It is not a total carbon-yield calculation when other monomers are consumed.

## Notebook execution

Install the packages and launch Jupyter as described in the [README](../README.md#getting-started). The kernel working directory should be `Scripts/`.

1. Start a fresh kernel and load the selected SBML file. Record the filename and repository revision.
2. Set the initial substrate bounds. The input `updated_5` model already contains `SK_bhb_c`; skip the cell that attempts to add this boundary again.
3. Choose whether to apply the later constraint cell before evaluating yields. Recalculate yields after changing bounds, since earlier stored values describe a different model state.
4. Use the inspection cells as needed. Skip the cell containing `cobra.io.write_sbml_model(model, sbml_filename)` unless intentionally exporting a model. Its default destination is the tracked `Models/260223_iJN1463_updated_6.sbml`; select a new filename for a new analysis.
5. Review the formula and polymer-reaction edits before applying them. They modify the in-memory model after the export cell.
6. Select the desired sweep implementation and verify its objective, composition, bounds, and output path. Confirm an optimal solver status before interpreting fluxes.

The saved execution counts are nonsequential, so existing cell outputs should be treated as historical results. A fresh sequential execution is not established as an exact reconstruction of the archived spreadsheet.

### Sweep implementation details

The final `set_polymer_uptake` helper only updates a sink when its calculated capacity is nonzero. When extending the analysis to mixtures with absent monomers, explicitly close those sinks to avoid retaining uptake bounds from a preceding scenario.

The final sweep's loop range uses the global `fixed_fraction_value`, while the function also accepts a `fixed_fraction` argument. Review both when changing the default half-mixture design.

### Output locations

Exports are relative to the kernel working directory. The two sweep cells write:

| Sweep cell | Output filename |
| --- | --- |
| Earlier, endpoint excluded | `250924_polymer_mixture_yield_sweep_half.xlsx` |
| Final, endpoint included | `251221_polymer_mixture_yield_sweep_half_inital_constrains.xlsx` |

The final cell's printed completion message names a different file. Use the path passed to `to_excel` to locate the output. Neither default export path points to the tracked spreadsheet in `Results/`.

## Archived data and results

### Single-monomer summary

[Data/250326_yield_Data.csv](../Data/250326_yield_Data.csv) contains four monomer-specific rows, each with a substrate uptake of 6 mmol gDW⁻¹ h⁻¹. It records growth, biomass yield, acetyl-CoA yield, R-3HB yield, and oxygen uptake. All four stored R-3HB yield values are zero.

Growth is labeled in h⁻¹ and biomass yield in gDW mmol⁻¹ substrate. The product-yield columns reuse biomass-yield units in their headers; their intended units require provenance verification before quantitative reuse. The current notebook does not read this CSV.

### Polymer-mixture summary

[Results/251218_polymer_mixture_yield_sweep_half.xlsx](../Results/251218_polymer_mixture_yield_sweep_half.xlsx) contains 153 scenarios, with columns identifying the fixed and varied polymers, mixture fractions, `mol/mol_yield`, and `product_flux`. Its scenario count matches the final sweep's default design. The exact generating model state is not recorded in a separate run manifest; matching the row count does not establish numerical reproduction. Account for the notebook's yield-definition issue when interpreting this historical export.

## Reporting an analysis

For an interpretable computational result, report the model filename and repository revision together with the Python, COBRApy, and solver versions; the objective and any biomass requirement; all modified reaction bounds; the feedstock composition and uptake capacity; the yield equation and units; and the output filename. Distinguish constraints stored in the SBML file from changes made during notebook execution.

Reaction-balance qualifications and differences among archived model snapshots are documented in the [model history](model-history.md).
