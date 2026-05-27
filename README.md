# BCR Lineage Tracer (B cell Linage Tree Maker)
In this project, I intend to create a Python tool for making B cell receptor (BCR) clonal lineage trees from single-cell sequencing data.
The tool infers the ancestral germline BCR, traces somatic hypermutation (SHM) events across cells, and outputs an annotated phylogenetic tree showing how B cell clones evolved and expanded.
B cell lineage tracing is a tool for understanding how the immune system builds and refines antibody responses. By capturing the order and location of these mutations, we can get insights into the selection shaping BCR evolution and antibody creation.

~The tool can be used in two ways: via a graphical user interface (GUI) for interactive exploration, or via the command line for scripted and batch workflows.~
Given paired heavy and light chain VDJ sequences, t

## Background: ##
* During an immune response, B cells in our body undergo clonal expansion and affinity maturation in order to produce a large pool of high‑affinity, antigen‑specific antibodies that can neutralize the pathogens more effectively.
As part of this process, the cells accumulate somatic hypermutations (SHM) in their BCR sequences, specifically in the V(D)J region, to improve the antigen.

* B cell lineage trees map the evolutionary "family history" of B cell clones as they mutate and divide during an immune response - all originated from one B cell.
  
* The trees are essential for couple of reasons:

  -**Tracking Affinity Maturation:** Trees illustrate how B cells repeatedly mutate and select for stronger, more precise antibody binding against specific pathogens.

  -**Isotype Switching:** By tracking the genetic changes, the trees can reveal class-switch recombination (change in the antiboy type), showing how B cells change their functional properties over time.

  -**Therapeutic Antibody Discovery:** We can use the trees to trace mutated B cell sequences backwards, identifying the most potent, broadly neutralizing antibodies to isolate for therapeutic use.

  -**Disease & Vaccine Research:** Analyzing tree shapes reveals whether an immune response is generating new, evolving antibodies or just re-stimulating older, less effective memory cells.
