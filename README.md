# Antivirulence Target Prioritization and Dual-Target Screening in *Klebsiella pneumoniae*

A computational drug-discovery project screening natural products against two virulence-associated and infection-relevant *Klebsiella pneumoniae* targets: **DsbA** (oxidative protein folding and virulence-factor maturation) and **Irp3** (yersiniabactin siderophore biosynthesis and iron acquisition).

The workflow combines subtractive genomics, structure prediction, structural validation, binding-pocket assessment, diversity-based natural-product selection, molecular docking, dual-target comparative analysis, ligand-efficiency prioritization, and ADMET/toxicity profiling.

---

## 1. Background and Target Identification

Targets were identified through subtractive genomics performed on the *K. pneumoniae* subsp. *pneumoniae* HS11286 reference proteome (**RefSeq assembly GCF\_000240185.1**), selected over the corresponding GenBank submission (GCA\_000240185.2) because the RefSeq assembly provided the curated reference annotation used for this analysis.

Proteome retrieval and initial redundancy reduction were performed upstream of the subtractive-genomics notebook. The raw RefSeq proteome containing 5,779 sequences was reduced using CD-HIT, producing 5,637 non-redundant clusters. The resulting `kp_nr.fasta` file was used as the starting input for the subtractive-genomics notebook.

From this point, the notebook performs the following workflow:

```text
kp_nr.fasta
    ↓
BLASTp against human proteome
    ↓
removal of human homologs
    ↓
BLASTp against reviewed UniProt/Swiss-Prot
    ↓
curated-reference filtering
    ↓
non-homologous candidate set
    ↓
annotation-based target prioritization
```

Swiss-Prot was used for the second filtering step instead of the originally planned Database of Essential Genes (DEG). The Swiss-Prot step should therefore **not** be interpreted as an essentiality filter; it removes proteins with strong matches in a curated and well-characterized reference set, which is a different criterion from predicted essentiality.

The remaining candidates were prioritized using annotation-based keywords including `dsba`, `siderophore`, `iron`, `hemolysin`, `secretion`, `porin`, and `transporter`.

This keyword-based step represents a **target-prioritization heuristic**, rather than experimental validation of essentiality or druggability.

### Selected targets

**DsbA (YP\_005224660.1)** is a thiol-disulfide oxidoreductase involved in oxidative protein folding and maturation of bacterial virulence-associated proteins. Its catalytic region contains the CXXC motif involving Cys110 and Cys113.

**Irp3 (YP\_005227768.1)** is associated with the yersiniabactin siderophore biosynthetic pathway and bacterial iron acquisition. The enzyme contains a catalytically relevant Tyr129 residue reported for this enzyme family.

DsbA and Irp3 are described here as **computationally prioritized candidate targets**, rather than newly discovered or experimentally validated drug targets.

At the time of the initial literature search, no clearly documented small-molecule inhibitor was identified for these specific *K. pneumoniae* targets. This search was not exhaustive and was not treated as proof that no such inhibitor exists.

A third candidate, **ClpV1**, was initially considered but was subsequently excluded because the extracted sequence was only 72 amino acids long, substantially shorter than expected for a full-length ClpV1 protein. DsbA and Irp3 therefore proceeded to structural analysis and compound screening.

---

## 2. Structure Prediction

Because suitable experimental structures were not used directly for the docking workflow, structures of DsbA and Irp3 were predicted using **ColabFold**, an AlphaFold2-based structure-prediction pipeline.

For DsbA, an experimental DsbA structure (PDB 4MCU) was considered as a potential template. However, sequence comparison produced weak and fragmented similarity and therefore did not provide a sufficiently reliable basis for direct template-based modeling.

ColabFold was consequently used for both targets.

| Target | Selected model    | pLDDT |   pTM |
| ------ | ----------------- | ----: | ----: |
| DsbA   | rank 001, model 2 |  91.4 | 0.799 |
| Irp3   | rank 001, model 5 |  94.5 | 0.921 |

Both models showed high predicted confidence. pLDDT values provide an indication of predicted local structural confidence, while pTM provides an estimate of global structural confidence.

### DsbA sequence length discrepancy

The ColabFold-predicted DsbA model was 261 aa, differing from the [CONFIRM: NCBI RefSeq / UniProt — which source?] annotation of 226 aa. SignalP-6.0 (Organism: Other, slow mode) identified a cleavable Sec/SPI signal peptide at the DsbA N-terminus (residues 1–24, cleavage probability 0.98), consistent with a manual von Heijne (-3,-1) rule check. Irp3 was predicted "OTHER" (no signal peptide) with 1.000000 confidence, consistent with its cytoplasmic role in yersiniabactin biosynthesis.

Since mature periplasmic DsbA does not retain its signal peptide, the signal peptide was trimmed prior to docking (mature DsbA: residues 25–261, 237 aa). The signal peptide was removed but remaining residues were not renumbered; residue numbers throughout this README therefore refer to the original precursor sequence numbering (e.g. Cys110/Cys113, Tyr112 below correspond directly to residue numbers in the prepared structure file).

SignalP prediction plots and result files are provided in `structure_validation/`.

The predicted structures and associated PAE/pLDDT outputs are provided in the `structures/` directory.

---

## 3. Structure Validation and Active-Site Assessment

The predicted models were assessed using **Ramachandran analysis** (SWISS-MODEL Structure Assessment) to evaluate backbone φ/ψ torsion-angle distributions. 

For DsbA, this was performed both before and after signal-peptide trimming:

| | Before trimming (261 aa) | After trimming (237 aa) |
|---|---|---|
| Ramachandran favored | 90.35% | **97.87%** |
| Ramachandran outliers | 8.11% | **0.85%** |
| MolProbity score | 2.35 | **1.34** |
| Clashscore | 7.16 | 5.63 |
| QMEAN4 z-score | -2.258 | **+1.419** |

26 of 28 bond-angle outliers in the untrimmed model localized to residues 1–30 — the signal peptide region — providing structural-quality evidence supporting the trimming decision, in addition to the biological rationale described in Section 2.

For Irp3 (confirmed signal-peptide-free by SignalP): Ramachandran favored 95.05%, outliers 1.65%, MolProbity score 2.14, clashscore 20.11 (scattered rather than localized; not re-assessed after PyRx energy minimization), QMEAN4 z-score +0.321. Verdict: passes, no trimming needed.

The corresponding validation plots are provided in:

```text
structure_validation/
├── dsba_ramachandran_before_trim.png
├── dsba_ramachandran_after_trim.png
├── irp3_ramachandran.png
├── dsba_signalp_plot.png
├── irp3_signalp_plot.png
├── dsba_signalp_results.txt
└── irp3_signalp_results.txt
```

### Protein Preparation

Discovery Studio Visualizer (free version) was used for hydrogen addition and geometry cleanup, since the paid-only "Prepare/Clean Protein" tools were unavailable. Atom counts were used to confirm Clean Geometry made only positional adjustments (no atom loss): DsbA 3,557 atoms (1,777 H), Irp3 5,720 atoms (2,841 H). Energy minimization was deferred to PyRx, performed automatically during receptor preparation.

Final prepared structures:
- `DsbA_prepped.pdb` (mature DsbA, residues 25–261)
- `irp3_prepped.pdb`

Binding pockets were subsequently evaluated using CASTpFold and cross-checked against functional residues, sequence information, and available structural literature.

### DsbA

The catalytic CXXC motif (Cys110-Val111-Tyr112-Cys113) was identified directly from sequence and confirmed present and structurally unaffected by trimming. The dominant CASTpFold pocket was approximately 723 Å² and included residues near this motif, including Tyr112. Full pocket-lining residue set: 62, 65, 66, 69, 70, 73, 112, 219–225, 242, 244–246.

The final docking grid was defined using the centroid of Cα atoms of selected pocket-lining residues.

### Irp3

Irp3 is a thiazolinyl imine reductase (Irp3/YbtU) in the yersiniabactin pathway. Meneely & Lamb (2016, PDB 5KVS) identify Tyr128 (*Y. enterocolitica* numbering) as the catalytic general acid. Pairwise alignment (Biopython) between the *K. pneumoniae* and *Y. enterocolitica* sequences revealed a +1 residue offset (extra N-terminal Met in the KP model), mapping the catalytic residue to Tyr129 in this structure — confirmed directly (residue 129 = TYR).

The Irp3 pocket was evaluated using available structural literature, sequence alignment, and CASTpFold analysis. CASTpFold identified the dominant pocket (~4,695 Å²) including Tyr129, alongside residues 17, 18, 21, 74–80, 101–103, 162–166, and others. Because this larger predicted surface also reflects NRPS-module-binding interactions rather than a compact small-molecule site, docking was focused on a smaller, functionally relevant region: residues 17, 18, 21, 101–103, 128–129, 162–166.

The final grid-box coordinates are provided in:

```text
results/docking/docking_grid_box_coordinates.json
```

---

## 4. Natural-Product Library Construction

Natural products were obtained from the **COCONUT database**.

The August 2026 release contained 738,827 natural-product compounds.

The library-processing workflow consisted of:

1. RDKit structure/SMILES validity checking
2. duplicate removal
3. biological-class filtering
4. Lipinski-based filtering
5. PAINS filtering
6. chemical-space characterization
7. diversity-based representative selection

After validity and duplicate filtering, 738,823 compounds remained.

A biological-class filter covering flavonoids, terpenoids, alkaloids, and phenolics was applied. Lipinski filtering was performed using a lenient criterion of no more than one violation to avoid excessive exclusion of natural products.

PAINS filtering removed 18,516 flagged compounds, resulting in **270,326 PAINS-clean compounds**.

---

## 5. Chemical-Space Diversity and Representative Selection

Docking the complete 270,326-compound library against two targets was computationally impractical.

Chemical-space diversity was therefore assessed using random 15,000-compound samples. The analysis indicated high structural diversity, with a mean pairwise Tanimoto distance of approximately 0.896 and substantial scaffold diversity.

The full library was subsequently reduced to **850 chemically diverse representative compounds** using the RDKit `MaxMinPicker` based on Tanimoto similarity.

This approach was selected instead of the initial k-means strategy because Euclidean distance is not an appropriate similarity measure for sparse binary molecular fingerprints.

PCA was used as a sanity check of the selected subset, while UMAP was used later for chemical-space visualization.

---

## 6. Docking Preparation and Benchmarking

Before large-scale docking, AutoDock Vina performance was benchmarked.

A benchmark involving 20 compounds against DsbA produced a mean docking time of approximately 63.3 seconds per compound at an exhaustiveness of 8.

Extrapolation to the complete 270,326-compound library against both targets indicated that full-library docking would require approximately 9,507 CPU hours, making it impractical within the available computational environment.

Therefore, the 850-compound diversity-representative subset was used for docking.

Ligand preparation successfully produced docking-ready structures for 834 of the 850 compounds.

---

## 7. Structure-Based Virtual Screening

AutoDock Vina was used to dock the prepared natural-product subset against both DsbA and Irp3.

The docking workflow was performed using the prepared receptor structures and target-specific grid boxes.

| Target | Successful docks | Excluded | Success rate |
| ------ | ---------------: | -------: | -----------: |
| DsbA   |        795 / 834 |       39 |        95.3% |
| Irp3   |        790 / 834 |       44 |        94.7% |

The excluded compounds were primarily associated with unsupported atom types or malformed PDBQT structures.

Individual docking results are provided in:

```text
results/docking/
├── docking_results_DsbA.csv
└── docking_results_Irp3.csv
```

---

## 8. Chemical-Space Score Overlay

UMAP using a Tanimoto/Jaccard-based metric was used to visualize the chemical space of the screened compounds and overlay docking scores.

Strong-scoring compounds were distributed across multiple regions of chemical space rather than being confined to a single dominant scaffold, suggesting chemical diversity among the prioritized compounds.

The resulting visualization is provided in:

```text
figures/chemical_space_score_overlay.png
```

The PCA analysis used as a preliminary sanity check is provided in:

```text
figures/pca_850_sanity_check.png
```

---

## 9. Dual-Target Comparative Analysis

A total of **790 compounds received valid docking scores against both targets**.

The raw DsbA and Irp3 docking scores showed a correlation of approximately **r = 0.887**.

However, molecular weight was also strongly correlated with the consensus docking score (**r = -0.746**), whereas rotatable-bond count showed a much weaker relationship (**r = -0.101**).

This diagnostic analysis indicated that **molecular size may have contributed substantially to the apparent agreement between the two docking scores**.

Therefore, raw docking-score correlation was not interpreted as evidence that the compounds genuinely possessed equivalent binding behavior toward both targets.

To reduce the influence of molecular size, ligand efficiency (LE) was calculated as:

```text
LE = -(docking affinity) / heavy-atom count
```

The mean ligand efficiency across both targets was used for comparative prioritization.

The prioritization procedure was:

1. Docking affinity ≤ −7.0 kcal/mol against both targets
2. Ranking by consensus ligand efficiency
3. Statistical filtering using the mean + 1 SD cutoff

This resulted in **19 compounds** for downstream ADMET analysis.

Detailed dual-target analysis is provided in:

```text
results/docking/dual_target_analysis_with_LE.csv
```

---

## 10. ADMET and Toxicity Screening

The 19 dual-target candidates were evaluated using:

- SwissADME
- pkCSM
- ProTox-3

Drug-likeness, pharmacokinetic characteristics, and toxicity-related endpoints were evaluated together.

Toxicity-related endpoints from pkCSM and ProTox-3 were evaluated alongside SwissADME drug-likeness and pharmacokinetic properties.

The baseline filtering criteria included:

- no predicted AMES mutagenicity
- no predicted hepatotoxicity
- Lipinski ≤1 violation
- acceptable predicted acute toxicity class

Following this filtering, **five ADMET-prioritized candidates** remained:

- CNP0085979.0
- CNP0329726.1
- CNP0383885.0
- CNP0167861.2
- CNP0374284.1

The compiled docking, ligand-efficiency, and ADMET results are provided in:

```text
results/admet/final_compiled_docking_LE_ADMET.csv
```

These candidates represent **computationally prioritized compounds**, not experimentally confirmed inhibitors.

---

## 11. Project Status and Limitations

**Project status:** Computational screening and prioritization was completed through docking, dual-target comparative analysis, ligand-efficiency assessment, and ADMET/toxicity profiling. Five compounds were retained as final computational candidates.

Molecular dynamics simulation was not performed as part of this repository.

The principal limitations are:

- Docking scores represent predicted binding affinity and are not experimental measurements.
- The DsbA and Irp3 structures used for docking were computationally predicted models.
- High pLDDT/pTM values indicate predicted structural confidence but do not substitute for experimental structural validation.
- Ramachandran analysis provides structural-quality information but does not experimentally validate the predicted structures.
- CASTpFold and literature/sequence-based pocket assessment provide computational support for the selected docking regions.
- The diversity analysis was based on representative samples rather than exhaustive pairwise comparison of the complete library.
- Apparent agreement between DsbA and Irp3 docking scores was found to be substantially influenced by molecular size.
- ADMET predictions are computational estimates and require experimental validation.
- No claim is made that the final compounds are experimentally confirmed DsbA/Irp3 inhibitors.

Molecular dynamics was considered a potential future extension, but it requires substantially greater computational resources and continuous trajectory generation, making it impractical within the available CPU-only environment.

---

## 12. Repository Structure

```text
DsbA-Irp3-virtual-screening/
│
├── data/
│   ├── kp_nr.fasta
│   ├── kp_non_human.fasta
│   ├── kp_final.fasta
│   └── kp_final_annotations.txt
│
├── figures/
│   ├── chemical_space_score_overlay.png
│   └── pca_850_sanity_check.png
│
├── notebooks/
│   ├── 01_subtractive_genomics.ipynb
│   └── 02_dual_target_screening.ipynb
│
├── results/
│   ├── admet/
│   │   └── final_compiled_docking_LE_ADMET.csv
│   │
│   ├── docking/
│   │   ├── docking_grid_box_coordinates.json
│   │   ├── docking_results_DsbA.csv
│   │   ├── docking_results_Irp3.csv
│   │   └── dual_target_analysis_with_LE.csv
│   │
│   └── sg/
│       ├── kp_vs_human.txt
│       ├── kp_vs_swiss.txt
│       ├── matched_ids.txt
│       ├── kp_priority_candidates.txt
│       ├── kp_tier1_candidates.txt
│       └── final_2_targets.fasta
│
├── structure_validation/
│   ├── dsba_ramachandran_before_trim.png
│   ├── dsba_ramachandran_after_trim.png
│   ├── irp3_ramachandran.png
│   ├── dsba_signalp_plot.png
│   ├── irp3_signalp_plot.png
│   ├── dsba_signalp_results.txt
│   └── irp3_signalp_results.txt
|
├── structures/
│   ├── dsba/
│   │   ├── colabfold/
│   │   └── prepared/
│   │
│   └── irp3/
│       ├── colabfold/
|       └── prepared/       
│
├── LICENSE
│
└── README.md
```

The `structures/` directory contains the predicted protein structures, ColabFold confidence outputs, and docking-prepared receptor structures.

The `structure_validation/` directory contains the Ramachandran validation results for the predicted DsbA and Irp3 structures.

---

## 13. Tools and Data Sources

### Sequence and target identification

- NCBI RefSeq
- BLASTp
- UniProt/Swiss-Prot
- CD-HIT

### Structure prediction and validation

- ColabFold / AlphaFold2
- Ramachandran analysis
- CASTpFold
- Biopython

### Compound processing and chemical-space analysis

- COCONUT Natural Products Database
- RDKit
- Open Babel
- scikit-learn
- UMAP

### Molecular docking

- AutoDock Vina

### ADMET and toxicity

- SwissADME
- pkCSM
- ProTox-3

### Computing

- Python
- Google Colab

External software and web-based resources may produce different results across versions or database releases. The versions available during the original analysis were not systematically recorded for every tool; therefore, this repository documents the workflow and outputs rather than claiming exact environment-level reproducibility.

---

## 14. License

The original code, notebooks, and analysis scripts in this repository are released under the **MIT License**.

Third-party databases, structures, software, and web-based prediction services remain subject to their respective licenses and terms of use.

The COCONUT natural-product database is provided under its stated CC0 terms.

---

## 15. Disclaimer

This repository presents a computational workflow for prioritizing potential natural-product candidates against two *K. pneumoniae* proteins.

The results should not be interpreted as evidence of experimentally confirmed inhibition, therapeutic efficacy, or clinical activity. Experimental biochemical and microbiological validation would be required to establish inhibitory activity and biological relevance.
