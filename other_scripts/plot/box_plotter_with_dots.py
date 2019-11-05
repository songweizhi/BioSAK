import statistics as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def box_plotter_with_dots(HGT_per_Mbp_list, HGT_num_per_genome_boxplot):

    HGT_per_Mbp_list_with_class = []
    for each in HGT_per_Mbp_list:
        HGT_per_Mbp_list_with_class.append([each, 1])

    HGT_per_Mbp_list_with_class_df = pd.DataFrame(np.array(HGT_per_Mbp_list_with_class), columns=['num', 'label'])

    sns.set(style="whitegrid")  # whitegrid

    # get boxplot
    orient = 'v'
    ax = sns.boxplot(x="label", y="num", data=HGT_per_Mbp_list_with_class_df, orient=orient, width=0.5,
                     showfliers=False)
    ax = sns.stripplot(x="label", y="num", data=HGT_per_Mbp_list_with_class_df, orient=orient, color='orange', size=3,
                       jitter=0.23)
    ax.set(ylabel='Number of HGT/Mbp sequence', xlabel='Genome')

    # plt.show()
    plt.savefig(HGT_num_per_genome_boxplot, dpi=300)
    plt.close()


box_plotter_with_dots([1,2,3,4,5], '/Users/songweizhi/Desktop/test.png')
