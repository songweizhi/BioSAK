import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector


# base = importr('base')
# utils = importr('utils')

utils = rpackages.importr('googleVis')


packages_needed = ['googleVis']
for each_package in packages_needed:
    if not rpackages.isinstalled(each_package):
        utils.install_packages(each_package)
    else:
        pass

# read in data
df = robjects.DataFrame.from_csvfile('/Users/songweizhi/Desktop/googleVis_input.txt')

sankey_plot = robjects.r['gvisSankey'](df, option = robjects.r['list'](sankey = "{node:{colorMode:'unique', labelPadding: 10 },link:{colorMode:'source'}}", height = 800, width = 600))

#robjects.r['plot'](sankey_plot)

out = open('/Users/songweizhi/Desktop/google.html', 'w')
out.write(str(sankey_plot))
out.close()
