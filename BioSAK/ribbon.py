import os
import glob
import gzip
import argparse
from Bio import SeqIO


ribbon_usage = '''
============================= ribbon example commands =============================

# Dependencies: opd

# example commands
odp="/scratch/PI/ocessongwz/Software/odp/scripts/odp"
ribbon="/scratch/PI/ocessongwz/Software/odp/scripts/odp_rbh_to_ribbon"
BioSAK ribbon -plot_lgs -co chrom_order.txt -o op_dir -fa demo_files -pep demo_files -chrom demo_files -f -odp $odp -odp_rbh_to_ribbon $ribbon

Note:
1. Species names can't have '_' char
2. Reference: https://github.com/conchoecia/odp
3. Demo data: https://github.com/songweizhi/BioSAK/tree/master/demo_data/ribbon

===================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def gff2chrom(gff_in, chrom_out):

    # This function was modified based on the NCBIgff2chrom.py from https://github.com/conchoecia/odp
    # This program parses a NCBI GFF annotation and generates a .chrom file
    # see https://github.com/conchoecia/odp for the specification

    gzipped = False
    for thisend in [".gz", ".gzip", ".GZ", ".GZIP", ".gzipped", ".GZIPPED"]:
        if gff_in.endswith(thisend):
            gzipped = True

    if gzipped:
        handle = gzip.open(gff_in, 'rt')
    else:
        handle = open(gff_in, "r")

    prots = dict()
    for line in handle:
        line = line.strip()
        splitd = line.split("\t")
        if line and len(splitd) > 7 and splitd[2] == "CDS" and "protein_id=" in line:
            pid = [x for x in splitd[8].split(";") if x.startswith("protein_id=")][0].replace("protein_id=", "")
            scaf = splitd[0]
            strand = splitd[6]
            start = int(splitd[3])
            stop = int(splitd[3])
            if pid not in prots:
                prots[pid] = {"scaf": scaf, "strand": strand, "start": start, "stop": stop}
            else:
                if start < prots[pid]["start"]:
                    prots[pid]["start"] = start
                if stop > prots[pid]["stop"]:
                    prots[pid]["stop"] = stop
    handle.close()

    # write out .chrom file
    chrom_out_handle = open(chrom_out, 'w')
    for pid in prots:
        chrom_out_handle.write("{}\t{}\t{}\t{}\t{}\n".format(pid, prots[pid]["scaf"], prots[pid]["strand"], prots[pid]["start"], prots[pid]["stop"]))
    chrom_out_handle.close()


def ribbon(args):

    fa_in                   = args['fa']
    fa_x                    = args['fax']
    gff_in                  = args['gff']
    gff_x                   = args['gffx']
    gbk_in                  = args['gbk']
    gbk_x                   = args['gbkx']
    pep_in                  = args['pep']
    pep_x                   = args['pepx']
    chrom_in                = args['chrom']
    chrom_x                 = args['chromx']
    op_dir                  = args['o']
    species_order_txt       = args['so']
    chrom_order_txt         = args['co']
    thread_num              = args['t']
    force_overwrite         = args['f']
    minscafsize             = args['m']
    odp_exe                 = args['odp']
    odp_rbh_to_ribbon_exe   = args['odp_rbh_to_ribbon']
    plot_all                = args['plot_all']
    plot_lgs                = args['plot_lgs']

    ################################################# define file name #################################################

    wd_path = os.path.abspath(os.getcwd())

    get_macrosynteny_plot_wd = '%s/get_macrosynteny_plot' % op_dir
    get_ribbon_diagram_wd    = '%s/get_ribbon_diagram'    % op_dir
    config_yaml_macrosynteny = '%s/config.yaml'           % get_macrosynteny_plot_wd
    cmd_txt                  = '%s/commands.txt'          % op_dir

    ####################################################################################################################

    fa_f_re   = '%s/*.%s' % (fa_in, fa_x)
    fa_f_list = glob.glob(fa_f_re)

    species_to_chrom_dict = dict()
    for each_fa in fa_f_list:
        _, _, fa_base, _ = sep_path_basename_ext(each_fa)
        qualified_ctg_set = []
        for each_seq in SeqIO.parse(each_fa, 'fasta'):
            seq_len = len(each_seq.seq)
            if seq_len >= minscafsize:
                qualified_ctg_set.append(each_seq.id)
        species_to_chrom_dict[fa_base] = qualified_ctg_set

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % get_macrosynteny_plot_wd)
    os.system('mkdir %s' % get_ribbon_diagram_wd)

    ####################################################################################################################

    species_order_list = []
    if species_order_txt is not None:
        for each_s in open(species_order_txt):
            species_order_list.append(each_s.strip())
    else:
        species_order_list = sorted(list(species_to_chrom_dict.keys()))

    chrom_order_dict = dict()
    if chrom_order_txt is not None:
        for each_c in open(chrom_order_txt):
            each_c_split = each_c.strip().split('\t')
            s_id = each_c_split[0]
            c_id = each_c_split[1]
            if s_id not in chrom_order_dict:
                chrom_order_dict[s_id] = []
            chrom_order_dict[s_id].append(c_id)

    fa_in_abspath                    = os.path.abspath(fa_in)
    op_dir_abspath                   = os.path.abspath(op_dir)
    get_macrosynteny_plot_wd_abspath = os.path.abspath(get_macrosynteny_plot_wd)
    get_ribbon_diagram_wd_abspath    = os.path.abspath(get_ribbon_diagram_wd)

    ################################### get config.yaml for getting macrosynteny plot ##################################

    # get config.yaml for getting macrosynteny plot
    config_yaml_macrosynteny_handle = open(config_yaml_macrosynteny, 'w')
    config_yaml_macrosynteny_handle.write('ignore_autobreaks: True\n')
    config_yaml_macrosynteny_handle.write('diamond_or_blastp: "diamond"\n')
    config_yaml_macrosynteny_handle.write('duplicate_proteins: "pass"\n')

    # write out plot_LGs
    if plot_lgs is True:
        config_yaml_macrosynteny_handle.write('plot_LGs: True\n')
    else:
        config_yaml_macrosynteny_handle.write('plot_LGs: False\n')

    config_yaml_macrosynteny_handle.write('plot_sp_sp: True\n')
    config_yaml_macrosynteny_handle.write('\n')

    # write out species section
    config_yaml_macrosynteny_handle.write('species:\n')
    for each_species in species_order_list:
        current_fa    = '%s/%s.fa'    % (fa_in_abspath, each_species)
        current_pep   = '%s/%s.pep'   % (fa_in_abspath, each_species)
        current_chrom = '%s/%s.chrom' % (fa_in_abspath, each_species)
        config_yaml_macrosynteny_handle.write('  %s:\n' % each_species)
        config_yaml_macrosynteny_handle.write('    proteins:    %s\n' % current_pep)
        config_yaml_macrosynteny_handle.write('    chrom:       %s\n' % current_chrom)
        config_yaml_macrosynteny_handle.write('    genome:      %s\n' % current_fa)
        config_yaml_macrosynteny_handle.write('    minscafsize: %s\n' % minscafsize)
    config_yaml_macrosynteny_handle.close()

    ############################################### get macrosynteny plot ##############################################

    snakemake_cmd_macrosynteny = 'snakemake --cores %s --snakefile %s' % (thread_num, odp_exe)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'w')
    cmd_txt_handle.write('cd %s\n' % get_macrosynteny_plot_wd_abspath)
    cmd_txt_handle.write(snakemake_cmd_macrosynteny + '\n')
    cmd_txt_handle.write('\n')
    cmd_txt_handle.close()

    # get macrosynteny plot
    print(snakemake_cmd_macrosynteny)
    os.chdir(get_macrosynteny_plot_wd_abspath)
    os.system(snakemake_cmd_macrosynteny)

    #################################### get config.yaml for getting ribbon diagram ####################################

    os.chdir(wd_path)
    step2_figures_abspath = '%s/odp/step2-figures'                % get_macrosynteny_plot_wd_abspath
    snakemake_cmd_ribbon  = 'snakemake --cores %s --snakefile %s' % (thread_num, odp_rbh_to_ribbon_exe)

    step2_figures_sub_dir_list = next(os.walk(step2_figures_abspath))[1]
    synteny_plot_list = ['synteny_nocolor']
    for each_sub_dir in step2_figures_sub_dir_list:
        if each_sub_dir.startswith('synteny_coloredby'):
            synteny_plot_list.append(each_sub_dir)

    for color_setting in synteny_plot_list:
        current_synteny_dir        = '%s/%s'          % (step2_figures_abspath, color_setting)
        current_wd                 = '%s/%s'          % (get_ribbon_diagram_wd_abspath, color_setting)
        current_config_yaml_ribbon = '%s/config.yaml' % current_wd
        os.mkdir(current_wd)

        with open(cmd_txt, 'a') as cmd_txt_handle:
            cmd_txt_handle.write('cd %s\n' % current_wd)
            cmd_txt_handle.write(snakemake_cmd_ribbon + '\n')
            cmd_txt_handle.write('\n')

        config_yaml_ribbon_handle = open(current_config_yaml_ribbon, 'w')
        config_yaml_ribbon_handle.write('chr_sort_order: optimal-chr-or\n')
        if plot_all is True:
            config_yaml_ribbon_handle.write('plot_all: True\n')
        else:
            config_yaml_ribbon_handle.write('plot_all: False\n')
        config_yaml_ribbon_handle.write('\n')

        # write out species_order section
        config_yaml_ribbon_handle.write('species_order:\n')
        for species in species_order_list:
            config_yaml_ribbon_handle.write('  - %s\n' % species)
        config_yaml_ribbon_handle.write('\n')

        # write out chromorder section
        config_yaml_ribbon_handle.write('chromorder:\n')
        if chrom_order_txt is not None:
            for each_s in chrom_order_dict:
                config_yaml_ribbon_handle.write('  %s:\n' % each_s)
                for each_chrom in chrom_order_dict[each_s]:
                    config_yaml_ribbon_handle.write('    - %s\n' % each_chrom)
            config_yaml_ribbon_handle.write('\n')
        else:
            for each_species in species_order_list:
                current_chrom_list = species_to_chrom_dict.get(each_species, [])
                if len(current_chrom_list) > 0:
                    config_yaml_ribbon_handle.write('  %s:\n' % each_species)
                    for each_chrom in current_chrom_list:
                        config_yaml_ribbon_handle.write('    - %s\n' % each_chrom)
                config_yaml_ribbon_handle.write('\n')

        # write out rbh_directory section
        config_yaml_ribbon_handle.write('rbh_directory: %s\n' % current_synteny_dir)
        config_yaml_ribbon_handle.write('\n')

        # write out species section
        config_yaml_ribbon_handle.write('species:\n')
        for each_species in species_order_list:
            current_fa    = '%s/%s.fa'    % (fa_in_abspath, each_species)
            current_pep   = '%s/%s.pep'   % (fa_in_abspath, each_species)
            current_chrom = '%s/%s.chrom' % (fa_in_abspath, each_species)
            config_yaml_ribbon_handle.write('  %s:\n' % each_species)
            config_yaml_ribbon_handle.write('    proteins:    %s\n' % current_pep)
            config_yaml_ribbon_handle.write('    chrom:       %s\n' % current_chrom)
            config_yaml_ribbon_handle.write('    genome:      %s\n' % current_fa)
            config_yaml_ribbon_handle.write('    minscafsize: %s\n' % minscafsize)
        config_yaml_ribbon_handle.close()

    ################################################## run snakemake ###################################################

    # get ribbon diagram
    for color_setting in synteny_plot_list:

        current_wd              = '%s/%s'           % (get_ribbon_diagram_wd_abspath, color_setting)
        pwd_plot_file           = '%s/output.pdf'   % current_wd
        pwd_plot_file_relocated = '%s/%s.pdf'       % (op_dir_abspath, color_setting.replace('synteny_', 'ribbon_'))

        os.chdir(current_wd)
        os.system(snakemake_cmd_ribbon)

        if os.path.isfile(pwd_plot_file) is True:
            os.system('mv %s %s' % (pwd_plot_file, pwd_plot_file_relocated))

    ################################################### final report ###################################################

    print('Ribbon diagram exported to %s' % op_dir)
    print('Done!')

    ####################################################################################################################


if __name__ == '__main__':

    ribbon_parser = argparse.ArgumentParser(usage=ribbon_usage)
    ribbon_parser.add_argument('-o',                    required=True,                                  help='output directory')
    ribbon_parser.add_argument('-fa',                   required=True,                                  help='fa file directory')
    ribbon_parser.add_argument('-fax',                  required=False, default='fa',                   help='fa file extension, default: fa')
    ribbon_parser.add_argument('-gff',                  required=False,                                 help='gff file directory')
    ribbon_parser.add_argument('-gffx',                 required=False, default='gff',                  help='gff file extension, default: gff')
    ribbon_parser.add_argument('-gbk',                  required=False,                                 help='gbk file directory')
    ribbon_parser.add_argument('-gbkx',                 required=False, default='gbk',                  help='gbk file extension, default: gbk')
    ribbon_parser.add_argument('-pep',                  required=False,                                 help='pep file directory')
    ribbon_parser.add_argument('-pepx',                 required=False, default='pep',                  help='pep file extension, default: pep')
    ribbon_parser.add_argument('-chrom',                required=False,                                 help='chrom file directory')
    ribbon_parser.add_argument('-chromx',               required=False, default='chrom',                help='chrom file extension, default: chrom')
    ribbon_parser.add_argument('-so',                   required=False, default=None,                   help='species order in the ribbon diagram')
    ribbon_parser.add_argument('-co',                   required=False, default=None,                   help='chromosome order in the ribbon diagram')
    ribbon_parser.add_argument('-m',                    required=False, type=int, default=10,           help='minscafsize, default: 10')
    ribbon_parser.add_argument('-t',                    required=False, type=int, default=1,            help='number of core, default: 1')
    ribbon_parser.add_argument('-f',                    required=False, action="store_true",            help='force overwrite')
    ribbon_parser.add_argument('-plot_all',             required=False, action="store_true",            help='plot_all')
    ribbon_parser.add_argument('-plot_lgs',             required=False, action="store_true",            help='plot_LGs')
    ribbon_parser.add_argument('-odp',                  required=False, default='odp',                  help='executable file odp, default: odp')
    ribbon_parser.add_argument('-odp_rbh_to_ribbon',    required=False, default='odp_rbh_to_ribbon',    help='executable file odp_rbh_to_ribbon, default: odp_rbh_to_ribbon')
    args = vars(ribbon_parser.parse_args())
    ribbon(args)
