# Model provenance and revision history

[Return to the repository overview](../README.md) · [Methods and reproducibility](methods.md)

The repository contains six dated SBML snapshots of the modified *Pseudomonas putida* KT2440 model, iJN1463. This history is based on comparisons of the deposited species, reactions, stoichiometry, flux bounds, and objectives. Dates are interpreted from the `YYMMDD` filename prefixes.

The latest deposited model is [260223_iJN1463_updated_6.sbml](../Models/260223_iJN1463_updated_6.sbml). The [yield calculation notebook](../Scripts/Hao_yield_calculation.ipynb) initially loads version 5 and includes an export cell targeting the version 6 filename. Specify the snapshot and any subsequent constraint changes when reporting a simulation.

## Deposited snapshots

Reaction counts include exchange, sink, and biomass reactions. Metabolite counts distinguish compartments. All six snapshots contain 1,462 SBML gene products and retain the model identifier `iJN1463`.

| Snapshot | Date | Reactions | Metabolites | Principal changes relative to the preceding snapshot |
| --- | --- | ---: | ---: | --- |
| [Version 1](../Models/241218_iJN1463_updated_1.sbml) | 2024-12-18 | 2,928 | 2,153 | Earliest deposited C2-platform snapshot; contains `ACPH` and has `ACS` constrained to zero. The unmodified starting model is not included for a complete baseline comparison. |
| [Version 2](../Models/250211_iJN1463_updated_2.sbml) | 2025-02-11 | 2,930 | 2,155 | Adds terephthalate and its dihydrodiol intermediate, with `TEREDEOXY` and `TEREDEHYD`. |
| [Version 3](../Models/250324_iJN1463_updated_3.sbml) | 2025-03-24 | 2,939 | 2,160 | Adds ethylene glycol and 1,4-butanediol assimilation reactions, `AACOA3HYD`, and associated boundary reactions. Corrects terephthalate/intermediate formulas and `TEREDEHYD` stoichiometry. Blocks `AACOAT`, `BDH`, and glucose exchange. |
| [Version 4](../Models/250325_iJN1463_updated_4.sbml) | 2025-03-25 | 2,942 | 2,162 | Adds adipate, adipoyl-CoA, `ADPCOAH`, `ADPCOAR`, and `SK_adpac_c`. Makes `3OXCOAT` reversible, blocks `ALCD2x`, and closes the terephthalate sink. |
| [Version 5](../Models/250709_iJN1463_updated_5.sbml) | 2025-07-09 | 2,943 | 2,162 | Adds the 3-hydroxybutyrate sink `SK_bhb_c` and makes it the objective. Adjusts `ICDHyr`, `PPC`, and `SUCOAS`; blocks `PC` and `THRS`; switches the open monomer sink from adipate to ethylene glycol. |
| [Version 6](../Models/260223_iJN1463_updated_6.sbml) | 2026-02-23 | 2,943 | 2,162 | Reopens `PC`, blocks `THRA`, and makes `THRS` reversible. No species or reactions are added relative to version 5. |

## Added pathway modules

Identifiers below use the COBRApy convention; the SBML serialization prefixes reaction identifiers with `R_` and species identifiers with `M_`.

| Module | Reaction identifiers | Relevant added metabolites |
| --- | --- | --- |
| Acetate phosphorylation | `ACPH` | Uses existing acetate, acetyl phosphate, ATP, and ADP metabolites. |
| Terephthalate assimilation | `TEREDEOXY`, `TEREDEHYD` | `terepa_c`, `12di12ditere_c` |
| Ethylene glycol oxidation | `ETHYGLYCHYD` | `etheglycol_c`, `etheglycol_e` |
| 1,4-Butanediol assimilation through 4-hydroxybutyrate | `14BUTADEH`, `4HYDBUDEH`, `4HYDBUSUCDEH` | `14btdl_c`, `4hbutald_c`, `h4but_c` |
| Acetoacetyl-CoA reduction | `AACOA3HYD` | Uses the existing `(R)`-3-hydroxybutyryl-CoA metabolite, `3hbcoa__R_c`. |
| Adipate activation and oxidation | `ADPCOAH`, `ADPCOAR` | `adpac_c`, `adpcoa_c` |

`TEREDEOXY` encodes the oxygen-dependent terephthalate dioxygenase step. Its stored reaction name incorrectly repeats the dehydrogenase label; use the identifier and stoichiometry to distinguish it from `TEREDEHYD`.

## Version 5 and version 6 constraints

Bounds are written as `(lower, upper)` in mmol gDW<sup>−1</sup> h<sup>−1</sup>. Signs refer to the reaction direction encoded in the model. A zero interval blocks flux; these reaction constraints alone do not establish a corresponding genetic knockout.

| Reaction | Version 5 | Version 6 | Interpretation |
| --- | --- | --- | --- |
| `ACS` | `(0, 0)` | `(0, 0)` | Acetyl-CoA synthetase blocked. |
| `ACPH` | `(-1000, 1000)` | `(-1000, 1000)` | Acetate phosphorylation remains reversible. |
| `SUCOAS` | `(-1000, 0)` | `(-1000, 0)` | Only the reverse of the stored succinyl-CoA synthetase reaction is allowed. |
| `ICDHyr` | `(0, 1000)` | `(0, 1000)` | Isocitrate dehydrogenase restricted to the forward direction. |
| `PPC` | `(0, 1000)` | `(0, 1000)` | Phosphoenolpyruvate carboxylase allowed forward. |
| `PC` | `(0, 0)` | `(0, 1000)` | Pyruvate carboxylase restored in version 6. |
| `THRA` | `(-999999, 999999)` | `(0, 0)` | Threonine aldolase blocked in version 6. |
| `THRS` | `(0, 0)` | `(-1000, 1000)` | Threonine synthase reopened in both directions in version 6. |
| `AACOAT`, `BDH`, `ALCD2x` | `(0, 0)` | `(0, 0)` | All three reactions remain blocked. |
| `EX_glc__D_e` | `(0, 0)` | `(0, 0)` | Glucose exchange closed. |
| `SK_etheglycol_c` | `(-10, 0)` | `(-10, 0)` | Cytosolic ethylene glycol supply allowed. |
| `SK_terepa_c`, `SK_14btdl_c`, `SK_adpac_c` | `(0, 0)` | `(0, 0)` | These monomer sinks are closed in the deposited snapshots. |
| `SK_bhb_c` | `(0, 1000)` | `(0, 1000)` | 3-Hydroxybutyrate removal and default optimization objective. |

Versions 1–4 maximize `BIOMASS_KT2440_WT3`; versions 5–6 maximize `SK_bhb_c`. Loading a snapshot therefore also loads a particular objective and medium configuration. The notebook changes constraints during analysis, so the deposited defaults do not describe every simulated condition. Monomer uptake through a cytosolic sink is a modeling boundary condition rather than an explicit transport mechanism.

## Stoichiometric and annotation status

All species in the six deposited files have both chemical formula and charge attributes. In version 6, the added pathway reactions listed above conserve the encoded elements and charge. The terephthalate dehydrogenase imbalance present in version 2 is resolved in version 3.

A direct stoichiometric audit of version 6, summing product minus reactant compositions with an absolute tolerance of `1e-6`, nevertheless identifies residuals in nine reactions already present in version 1:

| Reaction | Nonzero residuals |
| --- | --- |
| `FE3PYOVDDR` | Charge: −1 |
| `PQQFEP` | H: +1 |
| `MEPCT_1`, `HDH`, `PPRGL` | H: +1; charge: +1 |
| `APPAT`, `ARGDI`, `IDPh_1` | H: −1; charge: −1 |
| `REPHACCOAT` | H: −2; O: −1 |

Exchange, sink, demand, and biomass reactions are excluded from this balance statement because they represent system boundaries or aggregate biomass formation. The residuals refer to the encoded formulas and stoichiometry; they do not identify the appropriate biological correction by themselves. Accordingly, the deposited model should not be described as universally mass- and charge-balanced. This audit does not validate reaction directionality, thermodynamics, or predicted yields.
