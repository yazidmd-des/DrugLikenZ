# DrugLikenZ: An Open-Source Graphical User Interface for Batch Evaluation of Rule-of-Five Compliance
Mahdaoui, Y., & Zakkoumi, H. (2026) 
[![DOI](https://zenodo.org/badge/1150980540.svg)](https://doi.org/10.5281/zenodo.18500549)

## Overview

**DrugLikenZ** is a Python-based graphical user interface (GUI) designed for the rapid evaluation of chemical libraries against established drug-likeness criteria. The application automates property calculation and rule-of-five compliance checking, providing researchers with an intuitive visual heatmap of their results. By integrating local RDKit calculations with optional PubChem database queries, it offers a high-throughput workflow for preliminary drug discovery and chemical education.

## Technical Abstract

The core concept of this application is to rapidly predict the drug-likeness of chemical compounds by applying **Lipinski's Rule of Five (Ro5)**, a set of physicochemical criteria predictive of oral absorption and permeability. 

The calculation method involves using **RDKit** and querying the **PubChem database** via the compound SMILES string to retrieve or calculate key properties:
* **MW**: Molecular Weight
* **HBA**: Hydrogen Bond Acceptor count
* **HBD**: Hydrogen Bond Donor count
* **logP**: Lipophilicity

The application performs a binary evaluation (Pass/Fail) against specific thresholds. A compound is considered drug-like (accepted) if it violates **no more than one** of these four rules, meaning its total compliance score must be 3 or 4.

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

> Mahdaoui, Y., & Zakkoumi, H. (2026). **DrugLikenZ: An Open-Source Graphical User Interface for Batch Evaluation of Rule-of-Five Compliance**. Zenodo. ```https://doi.org/10.5281/zenodo.18500549```

[![DOI](https://zenodo.org/badge/1150980540.svg)](https://doi.org/10.5281/zenodo.18500549)

## License

This project is licensed under the **Creative Commons Attribution 4.0 International** 
```The Creative Commons Attribution license allows re-distribution and re-use of a licensed work on the condition that the creator is appropriately credited.```
