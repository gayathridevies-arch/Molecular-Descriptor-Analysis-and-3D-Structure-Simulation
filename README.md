# Molecular-Descriptor-Analysis-and-3D-Structure-Simulation
Python workflow for molecular descriptor calculation, drug-likeness filtering, property analysis, and 3D structure generation.
# Computational Molecular Property Analysis & 3D Structural Modeling Toolkit

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gayathridevies-arch/Molecular-Descriptor-Analysis-and-3D-Structure-Simulation
/blob/main/Molecular Descriptor Analysis and 3D Structure Simulation.ipynb)

A Python-based computational chemistry workflow that processes molecular blueprints (SMILES), calculates physicochemical descriptors, applies medicinal chemistry screening rules (Lipinski's Rule of Five and Veber's criteria), performs statistical data analysis, and executes 3D molecular geometry generation with interactive visualization.

---

## Project Overview

This repository demonstrates an end-to-end computational pipeline for evaluating small molecule drug-likeness and exploring molecular geometry. By bridging tabular chemical data with physical-chemical property evaluation and force-field coordinate generation, this script highlights core competencies in molecular modeling, data analysis, and software engineering.

### Key Scientific Capabilities
* **Structural Parsing & Conversion:** Converts 1D SMILES strings into rigorous RDKit graph representations and explicit 3D atomic coordinates.
* **Property Calculation:** Evaluates critical descriptors including Molecular Weight, Octanol-Water Partition Coefficient ($LogP$), Hydrogen Bond Donors/Acceptors, Topological Polar Surface Area ($TPSA$), Rotatable Bonds, and Fraction $sp^3$.
* **Compliance Filtering:** Automatically screens chemical datasets against established pharmaceutical guidelines (**Lipinski's Rule of Five** and **Veber's Rules**) to determine oral bioavailability potential.
* **Statistical Data Science:** Implements correlation matrices, distribution density histograms, and multivariable scatter plots to analyze property relationships across large chemical libraries.
* **3D Geometry & Visualization:** Generates 3D spatial conformations, executes energy minimization protocols via force fields, and renders interactive 3D molecular models using standardized PDB protocols.

---

## Technical Stack & Libraries

* **Programming Language:** Python 3.x
* **Molecular Processing:** RDKit (`Chem`, `AllChem`, `Descriptors`)
* **Data Manipulation & Analysis:** Pandas, NumPy
* **Statistical & Multi-panel Visualization:** Seaborn, Matplotlib (Object-oriented `ax` management)
* **Interactive 3D Graphics:** `py3Dmol` (WebGL-accelerated browser rendering via Google Colab)

---

## Pipeline Architecture

1. **Dataset Ingestion & Feature Extraction:** Loads molecular libraries (such as the Delaney ESOL dataset) and maps SMILES strings to computational objects to derive physicochemical descriptors.
2. **Rule-Based Filtering:** Applies logical masks to flag compounds exceeding molecular weight ($>500 \text{ g/mol}$), lipophilicity ($LogP > 5$), or hydrogen bonding thresholds.
3. **Exploratory Data Analysis (EDA):** Visualizes property distributions, density trends for drug-like vs. non-drug-like subspaces, and feature collinearity via correlation heatmaps.
4. **3D Conformation Generation:** Translates 2D topology into 3D space via distance geometry embedding (`EmbedMolecule`), adds explicit hydrogen atoms, and refines geometry using the Merck Molecular Force Field (`MMFFOptimizeMolecule`).
5. **Structural Export & Rendering:** Converts internal molecular coordinates into standard PDB text blocks for interactive inspection.

---

## Code Highlight: 3D Conformation Pipeline

```python
# Core 3D Modeling Pipeline Example for Aspirin
smiles_example = "CC(=O)OC1=CC=CC=C1C(=O)O"
mol = Chem.MolFromSmiles(smiles_example)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

# Export and interactive render bridge
pdb_block = Chem.MolToPDBBlock(mol)
view = py3Dmol.view(width=400, height=300)
view.addModel(pdb_block, "pdb")
view.setStyle({"stick": {}})
view.zoomTo()
view.show()


