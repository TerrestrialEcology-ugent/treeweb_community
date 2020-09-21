# Metadata

In this folder is gathered the different data files to reproduce the analyse presented in the manuscript. 

- plot_xy.csv, this file contains the coordinates of the 53 plots in Lambert72 (https://epsg.io/31370).
	- X: longitude
	- Y: latitude
	- id_plot: ID of the plot

- synthesis_expldata_raw.csv, the file contains information on the overstorey composition and the fragmentation of the plot. 
	- id_plot: ID of the plot
	- id_frag: ID of the forest fragment
	- fragm_std: standardized fragmentation index value, see Appendix B of the manuscript for more on the derivation of this index
	- speccomb: species composition of the overstorey, fsyl: Fagus sylvatica, qrob: Quercus robus, qrub: Quercus rubra, all: all 3 species together

- treeweb_*_format.csv, 9 files containing the community data. The first column is always the ID of the plot, the remaining columns are the species abundance. See Appendix C for the sampling procedure of the different communities.