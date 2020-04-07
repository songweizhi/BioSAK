import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


MAG_1 = ['MAG_1', 0.38, 'Kelp', 'Prot']
MAG_2 = ['MAG_2', 0.5, 'Tara', 'Prot']
MAG_3 = ['MAG_3', 0.6, 'Kelp', 'Gemm']
MAG_4 = ['MAG_4', 0.5, 'Kelp', 'Prot']
MAG_5 = ['MAG_5', 0.3, 'Tara', 'Gemm']
MAG_6 = ['MAG_6', 0.2, 'Kelp', 'Prot']
MAG_7 = ['MAG_7', 0.3, 'Tara', 'Prot']
MAG_9 = ['MAG_9', 0.38, 'Tara', 'Gemm']
MAG_10 = ['MAG_10', 0.2, 'Kelp', 'Gemm']
MAG_data = [MAG_1, MAG_2, MAG_3, MAG_4, MAG_5, MAG_6, MAG_7, MAG_9, MAG_10]

MAG_df = pd.DataFrame(MAG_data, columns=['MAG', 'Freq', 'Lifestyle', 'Phylum'])

box_order =   ['Kelp',   'Tara']
color_order = ['orange', 'lightblue']

sns.boxplot(data=MAG_df, x="Phylum", y="Freq",
            hue="Lifestyle", hue_order=box_order, palette=color_order)

plt.show()


