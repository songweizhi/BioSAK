# step_1, get dataframe with python script boxplot_matrix_COG.py
python3 boxplot_matrix_COG.py -in COG_func_stats -out COG_func_stats_df.txt -skip_1st_row -with_functional_description -in_percent

# step_2, important !!!
# please make sure sumarries for detected HGTs are in the last row in COG_func_stats_df.txt, as values in the last row will be used to position the coloured triangles and squares in the plot.

# step_3, get plot with Rscript COG_boxplot_last1row.R
Rscript COG_boxplot_last1row.R -i COG_func_stats_df.txt -o COG_func_stats_df.png
