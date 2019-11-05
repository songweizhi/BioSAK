import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# https://seaborn.pydata.org/generated/seaborn.swarmplot.html#seaborn.swarmplot


# # load example dataset
# sns.set(style="whitegrid")
# tips = sns.load_dataset("tips")
# print(tips)
#
# # get the boxplot
# ax = sns.boxplot(x="sex", y="total_bill", data=tips)
#
# # Use swarmplot() to show the datapoints on top of the boxes
# ax = sns.swarmplot(x="sex", y="total_bill", data=tips, color='blue')

# ax = sns.stripplot(x="label", y="num", data=HGT_per_Mbp_list_with_class_df, color='orange', size=3, orient='v', jitter=True)



# # show the plot
# plt.show()
#
#

#
#
# d = {'col1': [1, 2], 'col2': [3, 4]}
# df = pd.DataFrame(data=d)
# print(df)
#
# df2 = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), columns=['a', 'b', 'c'])
# print(df2)
#
#
df3 = pd.DataFrame(np.array([[1, 1], [3, 1]]), columns=['a', 'b'])
print(df3)
#
#
sns.set(style="whitegrid")
#
# get the boxplot
ax = sns.boxplot(x="b", y="a", data=df3)

# Use swarmplot() to show the datapoints on top of the boxes
ax = sns.swarmplot(x="b", y="a", data=df3, color='blue')

plt.show()

