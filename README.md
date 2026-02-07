
<div style="text-align: justify;">
  
# DrugLikenZ: An Open-Source Graphical User Interface for Multi-Rule Drug-Likeness and Batch Property Evaluation
  
Mahdaoui, Y., & Zakkoumi, H. (2026) <a href="https://doi.org/10.5281/zenodo.18500549"><img src="https://zenodo.org/badge/1150980540.svg" alt="DOI" style="vertical-align:middle;"></a>

<div align="center">
  <br>
  <a href="https://sourceforge.net/projects/druglikenz/files/latest/download">
    <img src="https://img.shields.io/badge/Download-DrugLikenZ-32CD32?style=for-the-badge&logo=sourceforge&logoColor=white" height="50" style="vertical-align:middle; margin-right:-5px;">
    <span style="background-color:#32CD32; display:inline-block; height:50px; width:50px; border-radius:0; vertical-align:middle; line-height:50px;">
      <img src="https://a.fsdn.com/allura/p/druglikenz/icon?8a6cbb5597b388c8f7e5e66210a6ab211fda8605ac88eb86a2cad663817f907e?&w=135" width="40" style="vertical-align:middle; margin-top:-5px;">
    </span>
  </a>
  <p><i>Version 1.0.2 | 91.9 MB via SourceForge</i></p>
</div>


## Overview

**DrugLikenZ** is a Python-based graphical user interface (GUI) designed for the rapid evaluation of chemical libraries against established **multi-rule drug-likeness criteria**. The application automates property calculation and compliance checking across multiple filters (Lipinski, Veber, Muegge, and RO3), providing researchers with an intuitive visual heatmap of their results. By integrating local **RDKit** calculations with optional **PubChem database** queries, it offers a high-throughput workflow for preliminary drug discovery and chemical education.

## Technical Abstract

> The core concept of this application is to rapidly predict the **drug-likeness** of chemical compounds by applying multiple established physicochemical filters, including **Lipinski's Rule of Five (Ro5)**, **Veber’s Rule**, **Muegge’s Method**, and the **Rule of Three (RO3)**. By evaluating parameters such as molecular weight, lipophilicity, flexibility, and surface area, the tool provides a comprehensive assessment of oral bioavailability, permeability, and lead-like potential for drug discovery.

### Calculation Method
The calculation method involves using **RDKit** and querying the **PubChem database** via the compound SMILES string to retrieve or calculate key molecular properties:

* **MW**: Molecular Weight
* **HBA**: Hydrogen Bond Acceptor count
* **HBD**: Hydrogen Bond Donor count
* **LogP**: Lipophilicity (Partition Coefficient)
* **ROTB**: Rotatable Bond count (Molecular Flexibility)
* **PSA**: Polar Surface Area
* **Rings**: Number of cyclic structures
* **Carbons**: Total Carbon atom count
* **Heteroatoms**: Count of non-carbon/non-hydrogen atoms
  
### Parameter Definitions
The drug-likeness of a molecule is evaluated using several physiochemical properties: **Molecular Weight (MW)** measures the total mass of atoms, influencing a drug's ability to cross biological membranes. **Hydrogen Bond Acceptors (HBA)** and **Donors (HBD)** count atoms capable of forming hydrogen bonds, which are critical for binding affinity and solubility. **LogP (Partition Coefficient)** determines lipophilicity, indicating the balance between water and fat solubility. **Rotatable Bonds (ROTB)** quantify molecular flexibility, where fewer bonds often correlate with better oral bioavailability. **Polar Surface Area (PSA)** calculates the surface sum of polar atoms to predict membrane permeability. **Rings**, **Carbons**, and **Heteroatoms** (non-C/H atoms like N or O) are used to quantify the structural complexity and chemical diversity of the compound.


The application performs a binary evaluation (Pass/Fail) against specific thresholds. A compound is considered drug-like (accepted) if it violates **no more than one** of these four rules, meaning its total compliance score must be 3 or 4.
  
 ### Rule Comparison Summary
This table compares the four screening methods available in **DrugLikenZ**. While Lipinski focuses on size and bonding, the Veber rule specifically targets flexibility (ROTB) and surface area (PSA).

| Parameter | Lipinski's Rule of Five | RO3 (Congreve et al.) | Muegge Method | Veber Rule |
| :--- | :--- | :--- | :--- | :--- |
| **MW** | <= 500 | < 300 | 200 to 600 | — |
| **LogP** | <= 5 | <= 3 | -2 to 5 | — |
| **HBA** | <= 10 | <= 3 | <= 10 | — |
| **HBD** | <= 5 | <= 3 | <= 5 | — |
| **ROTB** | — | <= 3 | <= 15 | <= 10 |
| **PSA** | — | <= 60 | <= 150 | <= 140 |
| **Rings** | — | — | <= 7 | — |
| **Carbons** | — | — | > 4 | — |
| **Heteroatoms** | — | — | > 1 | — |
| **Acceptance** | **Min. 3/4 pass** | **All 6/6 pass** | **All 9/9 pass** | **All 2/2 pass** |
  

## Key Features

* **Multiple Rulesets**: Supports Lipinski's Rule of Five, Miles Congreve et al. Rule of Three (RO3), Muegge Method, and Veber et al. rules.
* **Automated Data Retrieval**: Interfaces with the PubChem PUG REST API for compound metadata.
* **Visual Analytics**: Generates dynamic, color-coded heatmaps where green represents compliance and red represents violations.
* **Intelligent Formatting**: Names of compounds failing the selected criteria are automatically highlighted in red on the axis.
* **Large Dataset Management**: Includes a navigation system to handle files with hundreds of compounds.
* **Clean Export**: Save accepted candidates directly to a CSV file and export high-resolution heatmaps as PNG images.

## Requirements

The application requires Python 3.x and the following dependencies:
* `customtkinter`
* `rdkit`
* `pandas`
* `seaborn`
* `matplotlib`
* `pubchempy`
* `requests`

## Installation

1.  **Clone the repository**:
    ```bash
    git clone [https://github.com/yourusername/DrugLikenZ.git](https://github.com/yourusername/DrugLikenZ.git)
    ```
<br>
2.  **Install the necessary packages**:
    ```bash
    pip install customtkinter rdkit pandas seaborn matplotlib pubchempy requests
    ```

3.  **Run the application**:
    ```bash
    python main.py
    ```

## Usage

1.  Launch the application.
2.  Select your desired **filtering rule** from the dropdown menu.
3.  Use the **Browse** button to upload a CSV or TSV file containing SMILES strings.
4.  Specify the name of the **SMILES column** in the entry box.
5.  The application will automatically deduplicate and process the list.
6.  Use the navigation buttons to browse the compliance heatmaps.
7.  Click **Export** to save results or the heatmap image.

## Citation

If you use this software in your research, please use the following citation:

> Mahdaoui, Y., & Zakkoumi, H. (2026). **DrugLikenZ: An Open-Source Graphical User Interface for Multi-Rule Drug-Likeness and Batch Property Evaluation**. Zenodo. ```https://doi.org/10.5281/zenodo.18500549```
  
  **Check for updates here:** <a href="https://doi.org/10.5281/zenodo.18500549"><img src="https://zenodo.org/badge/1150980540.svg" alt="DOI" style="vertical-align:middle;"></a>

## License

This project is licensed under the **Creative Commons Attribution 4.0 International** 
```The Creative Commons Attribution license allows re-distribution and re-use of a licensed work on the condition that the creator is appropriately credited.```
  </div>