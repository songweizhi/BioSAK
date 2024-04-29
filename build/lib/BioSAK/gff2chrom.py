import gzip
import sys


def gff2chrom(gff_in, chrom_out):

    # This function was modified based on the NCBIgff2chrom.py from https://github.com/conchoecia/odp
    # This program parses a NCBI GFF annotation and generates a .chrom file
    # see https://github.com/conchoecia/odp for the specification

    gzipped = False
    for thisend in [".gz", ".gzip", ".GZ", ".GZIP", ".gzipped", ".GZIPPED"]:
        if gff_in.endswith(thisend):
            gzipped = True

    if gzipped:
        handle = gzip.open(gff_in,'rt')
    else:
        handle = open(gff_in, "r")

    prots = dict()
    for line in handle:
        line = line.strip()
        splitd=line.split("\t")
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


gff_in    = '/Users/songweizhi/Desktop/Limpet/GCF_932274485.2_xgPatVulg1.2_genomic.gff'
chrom_out = '/Users/songweizhi/Desktop/Limpet/GCF_932274485.2_xgPatVulg1.2_genomic.chrom'
gff2chrom(gff_in, chrom_out)

