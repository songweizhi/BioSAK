import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


MAG_df = pd.read_csv('Grouped_boxplot_with_extra_values.csv')

box_order =   ['Kelp',   'Planktonic']
color_order = ['orange', 'lightblue']

sns.boxplot(data=MAG_df, x="Phylum", y="Frequency",
            hue="Lifestyle", hue_order=box_order, palette=color_order)

sns.despine(offset=10, trim=True)

# add dots
plt.plot(0 - 0.2, 0.25, alpha=1, marker='v', markersize=10, markeredgewidth=0, color='deepskyblue')
plt.plot(0 + 0.2, 0.52, alpha=1, marker='^', markersize=10, markeredgewidth=0, color='coral')
plt.plot(1 - 0.2, 0.55, alpha=1, marker='^', markersize=10, markeredgewidth=0, color='coral')
plt.plot(1 + 0.2, 0.275, alpha=1, marker='v', markersize=10, markeredgewidth=0, color='deepskyblue')

plt.tight_layout()
plt.savefig('Grouped_boxplot_with_extra_values.png', bbox_inches='tight', dpi=300)
plt.close()


# # generate pandas dataframe and export to file
# MAG_1 = ['MAG_1', 0.38, 'Kelp', 'Proteobacteria']
# MAG_2 = ['MAG_2', 0.5, 'Planktonic', 'Proteobacteria']
# MAG_3 = ['MAG_3', 0.6, 'Kelp', 'Actinobacteriota']
# MAG_4 = ['MAG_4', 0.5, 'Kelp', 'Proteobacteria']
# MAG_5 = ['MAG_5', 0.3, 'Planktonic', 'Actinobacteriota']
# MAG_6 = ['MAG_6', 0.2, 'Kelp', 'Proteobacteria']
# MAG_7 = ['MAG_7', 0.3, 'Planktonic', 'Proteobacteria']
# MAG_9 = ['MAG_9', 0.38, 'Planktonic', 'Actinobacteriota']
# MAG_10 = ['MAG_10', 0.2, 'Kelp', 'Actinobacteriota']

# MAG_data = [MAG_1, MAG_2, MAG_3, MAG_4, MAG_5, MAG_6, MAG_7, MAG_9, MAG_10]
# MAG_df = pd.DataFrame(MAG_data, columns=['MAG', 'Frequency', 'Lifestyle', 'Phylum'])

# MAG_df.to_csv('/Users/songweizhi/Desktop/Grouped_boxplot_with_extra_values.csv', index=False)

