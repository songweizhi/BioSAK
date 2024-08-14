# step_1, get dataframe. Please note that the -with_functional_description option has been remove from boxplot_matrix_dbCAN.py
Python3 boxplot_matrix_dbCAN.py -in MAGs_dRep99_dbCAN_stats_GeneNumber -in_percent -skip_1st_row

# step_2, important !!!
# please make sure sumarries for detected HGTs are in the last row in COG_func_stats_df.txt, as values in the last row will be used to position the coloured triangles and squares in the plot.

# step_3, get plot
Rscript COG_boxplot_last1row.R -i MAGs_dRep99_dbCAN_stats_GeneNumber_CAZy_cate.txt -o MAGs_dRep99_dbCAN_stats_GeneNumber_CAZy_cate.png

