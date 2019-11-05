import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector


#base = importr('base')
#utils = importr('utils')

utils = rpackages.importr('utils')


packages_needed = ['hexbin']
for each_package in packages_needed:
    if not rpackages.isinstalled(each_package):
        utils.install_packages(each_package)
    else:
        pass


r_list = robjects.IntVector([1, 2, 3, 4, 5, 6, 7, 8])
r_matrix = robjects.r['matrix'](r_list, nrow = 2)

print(r_list)
print(r_matrix)


df = robjects.DataFrame.from_csvfile('/Users/weizhisong/Desktop/input_for_googlevis.csv')
print(df)

