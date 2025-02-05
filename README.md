This file contains the code to accompany the paper ``The Spine of a Supersingular $\ell$-isogeny Graph". 
Each Jupyter notebook filed (.ipynb) is a SageMath 10.4 file. 
Code to generate datasets was run on a MacBook Pro, Apple M3 Chip, 16 GB Memory, macOS Sonoma 14.6.1.
Runtimes are reported in the data generation files.

The file contents (in alphabetical order) are as follows:

* Center_DataGeneration.ipynb: File used to generate the data in center012925.csv
* Center_DataProcessing.ipynb: File used to generate graphs from the data in center012925.csv
* center012925.csv: Data for the 2-isogeny graph. Column header: 'p', 'size of spine', 'size of center','size of graph', 'Fp vertices in center'. The `graph' is the full Fp-bar 2-isogeny graph, and the center refers to the center of this graph.
* center_ell3_013125.csv: Similar to center012925.csv, but for 3-isogeny graph data instead of 2-isogeny graph data.
* diameter_data.csv: From SpineDiameters.ipynb, data on the diameters of the spine.
* graph_functions.sage: Functions for generating three supersingular elliptic curve $\ell$-isogeny graphs: The graph over $\mathbb{F}_p$, the spine, and the graph over $\overline{\mathbb{F}}_p$.
* Graph_Viz.ipynb: A file for visualizing the graphs generated in graph_functions.sage
* loops & multi-edges in G_l(Fp).pdf: A pdf listing pairs of primes $p,\ell$ for which the $\mathbb{F}_p$ supersingular $\ell$-isogeny graph has loops and/or multiedges
* Multi_Edges_in_G_ell_F_p.ipynb: Notebook used to compute data in the loops & multi-edges in G_l(Fp).pdf
* Small_Prime_Information.ipynb: Notebook to print out precise edge and vertex information for the $\mathbb{F}_p$ $\ell$-isogeny graph, and precisely how the edges and vertices change when mapped into the spine. 
* SmallCharacteristicGraphDescription.pdf: In this pdf, we explicitly compute the mapping from the $\mathbb{F}_p$ $2$- and $3$-isogeny graphs for certain small primes which are excluded from our general structure theorems.
* SpineDiameter_examples.ipynb: Graphs giving examples of the paper's theorems on spine diameters
* SpineDiameters.ipynb: Notebook for generating conglomerate data on the spine diameters. 
