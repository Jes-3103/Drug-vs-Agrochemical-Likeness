# Drug-vs-Agrochemical-Likeness
Python/RDKit screening of drugs vs. agrochemicals against Lipinski, Veber, Tice and PAINS filters

## Overview
Medicinal and agrochemical discovery both use rule-based filters to decide which compounds are worth screening. The best known is Lipinski's Rule of 5, derived from orally active pharmaceuticals and Tice´s rules later proposed an equivalent set for herbicides and insecticides.
This project applies both rule sets to the same 30 compounds (15 approved drugs and 15 commercial agrochemicals) and asks a simple question- do these filters actually distinguish the two classes?

## The pipeline:
1. Fetches canonical SMILES from PubChem by compound name
2. Calculates MW, cLogP, H-bond donors and acceptors, TPSA and rotatable bonds with RDKit
3. Applies Lipinski's Rule of 5, Veber's rules, Tice's agrochemical criteria and the PAINS substructure filter
4. Produces a results table and figures comparing the two classes

## Results
1. A plot of MW vs logP shows that agrochemicals occupy a smaller, less polar region than drugs. Mean molecular weight between the two categories differ by about 80 g/mol and mean TPSA by about 11 Å², while lipophilicity is essentially identical. Agrochemicals must cross waxy plant and insect cuticles, which favours smaller, less polar molecules, whereas orally administered drugs can tolerate more size and polarity. So this is consistent with their different uptake routes.

2. The agrochemicals pass the drug rule more often than the drugs do. This is because every Lipinski criterion is an upper limit, so small, simple molecules pass easily. On the other hand, Tice's criteria are harder for both classes because they are two-sided. Several of Tice's criteria are ranges, so a compound can fail for being too small or too polar as well as too large. Both classes therefore score lower against Tice than against Lipinski. Neither rule set works as a classifier. These rules describe where known actives have historically been found; they were never designed to discriminate between classes. PAINS filters are designed for screening libraries rather than marketed products; applying them here was a check on the method

3. Two of the fifteen drugs fail Lipinski's Rule of 5. Atorvastatin (MW 559, cLogP 6.3) breaks both the weight and lipophilicity limits, yet is one of the best-selling pharmaceuticals ever synthesised. Cyclosporin fails on size and hydrogen-bond count at MW 1203, and is a clinically essential immunosuppressant. In the other direction, lambda-cyhalothrin passes despite a cLogP of 6.5, because its molecular weight (450) and absence of hydrogen-bond donors keep it to a single violation — and one violation is tolerated. The rule is a statistical observation about orally absorbed compounds, not a criterion any individual molecule must satisfy.

## How to run
Create the environment and install the dependencies:
conda env create -f environment.yml
conda activate chemscreen

Then run the pipeline in order:
python fetch.py      # data/compounds.csv            -> data/compounds_with_smiles.csv
python analyse.py    # data/compounds_with_smiles.csv -> results/screening_results.csv
python figures.py    # results/screening_results.csv  -> results/*.png

fetch.py requires an internet connection for the PubChem lookups. Everything after that runs offline.

## Limitations
cLogP is calculated, not measured. RDKit uses the Wildman–Crippen atom-contribution method. Several compounds in the set (glyphosate, glufosinate, dicamba) are ionisable, so their behaviour at physiological or soil pH is better described by logD than logP.
Tice's criteria were defined using MLogP (Moriguchi), not Crippen cLogP. The two methods can differ by more than a log unit, so compounds near the 3.5 boundary may be classified differently here than in the original paper.
H-bond counts use RDKit's modern pattern-based definitions. Lipinski's original paper simply counted N + O atoms and NH + OH groups (Lipinski.NOCount and Lipinski.NHOHCount), which gives different values for some compounds.
TPSA counts only N and O by default, so the phosphorus in glyphosate and glufosinate contributes nothing. CalcTPSA(mol, includeSandP=True) includes it.
n = 30. This is a proof of concept. Statistically meaningful comparisons would require a few hundred to a few thousand compounds from ChEMBL.


## Environment notes
Building a working environment on Windows proved to be a substantial part of this project, and the diagnosis is worth recording:
1. Creating a figure (plt.subplots()) killed the Jupyter kernel with no Python error, while matplotlib imported cleanly and the same code ran without issue from the terminal. The fault was isolated by running matplotlib outside the notebook, which succeeded.
2. RDKit's drawing module failed separately with ImportError: DLL load failed while importing rdMolDraw2D, despite rdkit and librdkit being matched conda-forge builds (2026.09.2). This indicates a missing system-level dependency rather than a version mismatch.
3. As a result, figure generation and structure rendering currently run in separate conda environments. Resolving this to a single reproducible environment is outstanding.

## References
1. Lipinski, C. A. et al. Adv. Drug Deliv. Rev. 1997, 23, 3–25.
2. Veber, D. F. et al. J. Med. Chem. 2002, 45, 2615–2623.
3. Tice, C. M. Pest Manag. Sci. 2001, 57, 3–16.
4. Baell, J. B.; Holloway, G. A. J. Med. Chem. 2010, 53, 2719–2740. (PAINS)
