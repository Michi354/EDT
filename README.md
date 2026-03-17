# Evolutionary digital twin framework for optimal aminoglycoside dosing in neonates with suspected sepsis

## Main Entry Point

- `EcoMIC_v18_withValidation/runProject.m`

## Framework Overview

The project combines a mechanistic and data-driven workflow for neonatal
amikacin analysis and optimization. At a high level, `runProject.m`:

1. builds a virtual neonatal population with physiology-aware covariates,
2. evaluates the credibility of the generated population,
3. loads one of the clinical cohorts through a GUI-based cohort selector,
4. fits or loads serum creatinine prediction models,
5. derives patient-level physiology such as height, eGFR, cardiac output,
   and amikacin distribution volumes,
6. injects those quantities into the PBPK/PD model,
7. evaluates PK accuracy against observed amikacin concentrations,
8. compares the creatinine framework against the Krzyzanski reference model,
9. calibrates eco-evolutionary PD parameters on virtual and real cohorts,
10. optimizes amikacin regimens under PK/PD constraints.

The codebase therefore spans virtual population generation, mixed-effects
creatinine modeling, physiology updates, SimBiology simulation, bacterial PD
calibration, and regimen optimization.

## Repository Structure

- `EcoMIC_v18_withValidation/runProject.m`
  Main orchestrator for the complete workflow.
- `EcoMIC_v18_withValidation/PBPK_PD_EcoEvolutive.sbproj`
  SimBiology PBPK/PD project loaded and updated by the MATLAB scripts.
- `EcoMIC_v18_withValidation/dataverse_files/DATA_ACCESS.txt`
  Public data-access note replacing the distributed cohort spreadsheets.
- `EcoMIC_v18_withValidation/SavedModels/MODEL_REGENERATION.txt`
  Note describing which mixed-effects model files must be regenerated locally
  from approved data.

## PBPK/PD Project

The file `PBPK_PD_EcoEvolutive.sbproj` is the central SimBiology project used
by the framework. The MATLAB scripts update this model on a patient-by-patient
basis using cohort-specific or synthetic physiology.

From the code path in `runProject.m`, the PBPK/PD project is used to:
- load the base mechanistic amikacin model,
- update compartment capacities, partition coefficients, blood flows, and
  clearance from predicted physiology,
- apply time-varying physiology events when needed,
- simulate serum amikacin exposure,
- track eco-evolutionary bacterial states such as sensitive bacteria `Ns`,
  resistant bacteria `Nr`, and the adaptive state `u`,
- support downstream PD calibration and schedule optimization.

In practice, the repository couples learned creatinine models with a
mechanistic PBPK layer and an eco-evolutionary PD layer, rather than treating
PK and PD as disconnected steps.

## Clinical Datasets

The code uses three neonatal cohorts as primary data sources:
- `NONMEM.xlsx`
- `dataset as published creacool.xlsx`
- `data as published Amicool.xlsx`

Clinical characteristics of the three datasets used are summarized in Table 7,
"Clinical characteristics of patients in the three neonatal cohorts," in the
paper.

## Data Availability

The public repository does not include the patient-level cohort spreadsheets.


## Requirements

To run the workflow, the project expects a MATLAB R2024 environment with support for
the components exercised by the scripts, including:
- SimBiology v23.2
- Statistics and Machine Learning functionality
- Optimization functionality

The workflow also relies on interactive dataset selection through
`selectPipelineFromGUI.m`.

## Running The Workflow

1. Open MATLAB in `EcoMIC_v18_withValidation`.
2. Restore the required cohort spreadsheets locally if you have approved data
   access.
3. Re-train neural network models from the approved cohort data.
6. Run `runProject`.
