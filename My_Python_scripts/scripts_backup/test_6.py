import os
import glob
import shutil
from Bio import SeqIO
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages



def GvisSankey_plotter(pwd_googlevis_input, pwd_plot_html) :
    out = open(pwd_plot_html, 'w')
    utils = rpackages.importr('googleVis')
    packages_needed = ['googleVis']
    for each_package in packages_needed :
        if not rpackages.isinstalled(each_package) :
            utils.install_packages(each_package)
        else :
            pass

    df = robjects.DataFrame.from_csvfile(pwd_googlevis_input)
    sankey_plot = robjects.r['gvisSankey'](df,
                                           option = robjects.r['list'](
                                               sankey = "{node : {colorMode: 'unique', labelPadding: 10}, link:{colorMode: 'source'}}",
                                               height = 1800,
                                               width = 600))
    out.write(str(sankey_plot))
    out.close()


pwd_googlevis_input = '/Users/weizhisong/Desktop/222/output_MetaBATsuperspecific_vs_MyCC56mer/input_for_googlevis.csv'
pwd_plot_html = '/Users/weizhisong/Desktop/222/output_MetaBATsuperspecific_vs_MyCC56mer/googleVis_MetaBATsuperspecific_against_MyCC56mer.html'

GvisSankey_plotter(pwd_googlevis_input, pwd_plot_html)