import os
import argparse


RunGraphMB_usage = '''
=============== RunGraphMB example commands ===============

BioSAK RunGraphMB -gfa assembly_graph.gfa -o output_dir

===========================================================
'''


def gfa_to_depth(gfa_file, edge_fasta, depth_txt, depth_txt_MetaBAT):

    # parse gfa file
    edge_id_set = set()
    edge_seq_dict = dict()
    edge_len_dict = dict()
    edge_depth_dict = dict()
    for each_line in open(gfa_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('S'):
            edge_id = each_line_split[1]
            edge_seq = each_line_split[2]
            edge_len = len(edge_seq)

            try:
                edge_depth = int(each_line_split[3].split(':')[-1])
            except:
                edge_depth = float(each_line_split[3].split(':')[-1])

            edge_id_set.add(edge_id)
            edge_len_dict[edge_id] = edge_len
            edge_seq_dict[edge_id] = edge_seq
            edge_depth_dict[edge_id] = edge_depth

    # write out edge sequence
    edge_fasta_handle = open(edge_fasta, 'w')
    for each_edge in sorted([i for i in edge_id_set]):
        edge_seq = edge_seq_dict[each_edge]
        edge_fasta_handle.write('>%s\n' % each_edge)
        edge_fasta_handle.write('%s\n' % edge_seq)
    edge_fasta_handle.close()

    # write out depth info
    depth_txt_handle = open(depth_txt, 'w')
    depth_txt_handle.write('Sequence\tLength(bp)\tDepth(x)\n')
    for each_edge in sorted([i for i in edge_id_set]):
        edge_len = edge_len_dict[each_edge]
        edge_depth = edge_depth_dict[each_edge]
        depth_txt_handle.write('%s\t%s\t%s\n' % (each_edge, edge_len, edge_depth))
    depth_txt_handle.close()

    # write out depth info (MetaBAT format)
    depth_txt_MetaBAT_handle = open(depth_txt_MetaBAT, 'w')
    depth_txt_MetaBAT_handle.write('contigName\tcontigLen\ttotalAvgDepth\tsample1\tsample1-var\n')
    for each_edge in sorted([i for i in edge_id_set]):
        edge_len = edge_len_dict[each_edge]
        edge_depth = edge_depth_dict[each_edge]
        fake_depth_var = float("{0:.2f}".format(edge_depth*0.01))
        depth_txt_MetaBAT_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_edge, edge_len, edge_depth, edge_depth, fake_depth_var))
    depth_txt_MetaBAT_handle.close()


def RunGraphMB(arg_dict):

    gfa_file = arg_dict['gfa']
    op_dir   = arg_dict['o']

    # create folder
    if os.path.isdir(op_dir):
        os.system('rm -r %s' % op_dir)
    os.system('mkdir %s' % op_dir)

    # 1. copy gfa_file into GraphMB_input_dir
    os.system('cp %s %s/assembly_graph.gfa' % (gfa_file, op_dir))

    # 2. prepare depth file and edge sequence file
    depth_file_tmp = '%s/assembly_depth.tmp.txt'    % op_dir
    depth_file     = '%s/assembly_depth.txt'        % op_dir
    sequence_file  = '%s/assembly.fasta'            % op_dir
    gfa_to_depth(gfa_file, sequence_file, depth_file_tmp, depth_file)
    #os.system('rm %s' % depth_file_tmp)

    # 3. command for running GraphMB
    print('# Example commands for running GraphMB')
    cmd_GraphMB = 'conda activate graphmb\ngraphmb --assembly %s --outdir GraphMB_output_dir' % op_dir
    print(cmd_GraphMB)


if __name__ == '__main__':

    RunGraphMB_parser = argparse.ArgumentParser(usage=RunGraphMB_usage)
    RunGraphMB_parser.add_argument('-gfa', required=True, help='gfa file')
    RunGraphMB_parser.add_argument('-o',   required=True, help='output folder (i.e., input folder to GraphMB)')
    args = vars(RunGraphMB_parser.parse_args())
    RunGraphMB(args)
