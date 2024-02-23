source activate cellphonedb
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --iterations=1000 --threads=40 --counts-data hgnc_symbol
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_meta.txt 

