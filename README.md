Scripts to produce outputs for the microbial ocean atlas (see document `Project_Background.md` in this repository for more details)

Scripts with the prefix `SSI_` are used internally at CSRIO to calculate indices using the Australian microbiome primary data tables. Output data tables from the workflow are contained in the data directory in this repository. An overview of the workflow to generate the data products is described in the document ` script_workflow.md` in this repository.

Taxa, Traits or genes of interest can be included in the analysis by incorporating them into the `taxa_request.list`, ` trait_request.list` and ` gene_request.csv` documents respectively.

R scripts for the generation of species and community indices for a variety of environmental parameters are used as described in:
Brown, M.V., Ostrowski, M., Messer, L.F. et al. A marine heatwave drives significant shifts in pelagic microbiology. Commun Biol 7, 125 (2024). https://doi.org/10.1038/s42003-023-05702-4

Details on producing indices from ASV tables generated/downloaded through the Australian Microbiome data portal (https://data.bioplatforms.com/bpa/otu) are under development.

