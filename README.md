# BCR Lineage Tracer (B cell Linage Tree Maker)
In this project, I intend to create a Python tool for making B cell receptor (BCR) clonal lineage trees from single-cell sequencing data.
The tool infers the ancestral germline BCR, traces somatic hypermutation (SHM) events across cells, and outputs an annotated phylogenetic tree showing how B cell clones evolved and expanded.
B cell lineage tracing is a tool for understanding how the immune system builds and refines antibody responses. By capturing the order and location of these mutations, we can get insights into the selection shaping BCR evolution and antibody creation.


## Background: ##
* During an immune response, B cells in our body undergo clonal expansion and affinity maturation in order to produce a large pool of high‑affinity, antigen‑specific antibodies that can neutralize the pathogens more effectively.
As part of this process, the cells accumulate somatic hypermutations (SHM) in their BCR sequences, specifically in the V(D)J region, to improve the antigen.

* B cell lineage trees map the evolutionary "family history" of B cell clones as they mutate and divide during an immune response - all originated from one B cell.
  
* The trees are essential for couple of reasons:

  -**Tracking Affinity Maturation:** Trees illustrate how B cells repeatedly mutate and select for stronger, more precise antibody binding against specific pathogens.

  -**Isotype Switching:** By tracking the genetic changes, the trees can reveal class-switch recombination (change in the antiboy type), showing how B cells change their functional properties over time.

  -**Therapeutic Antibody Discovery:** We can use the trees to trace mutated B cell sequences backwards, identifying the most potent, broadly neutralizing antibodies to isolate for therapeutic use.

  -**Disease & Vaccine Research:** Analyzing tree shapes reveals whether an immune response is generating new, evolving antibodies or just re-stimulating older, less effective memory cells.


## The tool: ##
### How is it working:
I developed this Python script to construct BCR clonal lineage trees based on the fundamental principle of clonal evolution. The tool maps the evolutionary "family history" of B cell clones by rooting each tree in its unmutated common ancestor (germline)—representing the cell’s original state before it encountered an antigen. As these B cells divide during an immune response, they accumulate unique genetic changes through the SHM process. By analyzing the hierarchy of shared mutations, the script reconstructs the branching order of the cells, where the length of each branch reflects the number of genetic changes acquired over time. To ensure high-fidelity results, the code incorporates critical biological constraints such as isotype switching (the irreversible transition between antibody types) and tissue connectivity, allowing for a professional visualization of how immune lineages mature and spread throughout different organs.

### Input:
The tool will take an **input of tabular data in .xlsx /.csv /.tsv**  with one row per cell.
The required columns will be:

| Column | Description |
|---|---|
| `cell_id` | Unique cell barcode |
| `clone_id` | Clone group identifier |
| `clone_count` | Number of cells in this clone |
| `VDJ_sequence_H` | Heavy chain VDJ nucleotide sequence |
| `VDJ_sequence_L` | Light chain VDJ nucleotide sequence |
| `VDJ_aa_sequence_H` | Heavy chain VDJ amino acid sequence |
| `VDJ_aa_sequence_L` | Light chain VDJ amino acid sequence |
| `v_call_h` | Heavy chain V gene call |
| `d_call_h` | Heavy chain D gene call |
| `j_call_h` | Heavy chain J gene call |
| `v_call_l` | Light chain V gene call |
| `j_call_l` | Light chain J gene call |
| `mu_count_h` | Somatic hypermutation count, heavy chain |
| `mu_count_l` | Somatic hypermutation count, light chain |
| `cluster_annotated` | Cell type annotation |
| `c_call` | Isotype |
| `sample_id` | Sample / timepoint identifier |

* I decided to start with processed files rather than raw data because raw scBCR‑seq output typically comes in complex bioinformatics formats, and handling it can easily become a full project on its own. For this project, I prefer to focus on a less explored biological insights rather than reconstructing known piplens that are very common.

### Output:
**The final output will be a structured folder with the image of the tree and table files.**

| Output | Format | Description |
|---|---|---|
| Lineage tree | PNG/JPG | Phylogenetic tree with nodes colored by cell type or isotype and edges or notes to sign mutations |
| Mutation table | xlsx/CSV | List of nucleotide and amino acid changes |

Example of the tree structure:

<img width="425" height="650" alt="image" src="https://github.com/user-attachments/assets/752dfbb7-862a-4431-87c5-cc24752ede08" />



## The technicalities: ##
As for now, the requirements.txt will include:
* pandas
* openpyxl
* biopython
* ete3
* matplotlib
* scipy
* numpy
* tkinter (for GUI )

## Installation and running the project

**Requirements:** Python 3.9+

### 1. Install dependencies

pip install -r requirements.txt

### 2. Run the project

**Option A — Graphical Interface (recommended):**

python gui.py

This opens a window where you can upload your `.xlsx` / `.csv` / `.tsv` file,
select the clone to analyze, and view or export the results.

**Option B — Command Line:**

python run_lineage.py - input your_file.xlsx

### 3. Run the tests

pytest tests/


### Notes:
This project is part of the Python Programming Course at the Weizmann Institute of Science.

You can view the course main repository here: https://github.com/Code-Maven/wis-python-course-2026-03






