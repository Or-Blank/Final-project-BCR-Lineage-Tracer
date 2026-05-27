# BCR Lineage Tracer (B cell Linage Tree Maker)
In this project, I intend to create a Python tool for making B cell receptor (BCR) clonal lineage trees from single-cell sequencing data.
The tool infers the ancestral germline BCR, traces somatic hypermutation (SHM) events across cells, and outputs an annotated phylogenetic tree showing how B cell clones evolved and expanded.
B cell lineage tracing is a tool for understanding how the immune system builds and refines antibody responses. By capturing the order and location of these mutations, we can get insights into the selection shaping BCR evolution and antibody creation.

~The tool can be used in two ways: via a graphical user interface (GUI) for interactive exploration, or via the command line for scripted and batch workflows.~
Given paired heavy and light chain VDJ sequences, t

# Background: #
During an immune response, B cells in our body undergo clonal expansion and affinity maturation in order to produce a large pool of high‑affinity, antigen‑specific antibodies that can neutralize the pathogens more effectively.
As part of this process, the cells cumulate somatic hypermutations (SHM) in their BCR sequences to improve the antigen.

Single-cell BCR sequencing captures a snapshot of this process, but the evolutionary relationships between cells are not immediately visible in the raw data.

B cell lineage trees map the evolutionary "family history" of B cell clones as they mutate and divide during an immune response. They are essential for tracking how effectively our bodies fight infections, developing targeted vaccines, and understanding what goes wrong in autoimmune diseases and blood cancers

Tracking Affinity Maturation: Trees illustrate how B cells repeatedly mutate and select for stronger, more precise antibody binding against specific pathogens.
Isotype Switching: By tracking genetic changes, trees reveal class-switch recombination (e.g., transitioning from IgM to IgG), showing how B cells change their functional properties over time.
Therapeutic Antibody Discovery: Scientists use these trees to trace mutated B cell sequences backward, identifying the most potent, broadly neutralizing antibodies to isolate for therapeutic use.
Disease & Vaccine Research: Analyzing tree shapes reveals whether an immune response is generating new, evolving antibodies or just re-stimulating older, less effective memory cells
