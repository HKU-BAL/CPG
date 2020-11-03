import os
from Bio import SeqIO
import argparse

def pass_contigs(Ensure_contigs, pass_contigs):
    F_contigs = open(Ensure_contigs, 'r')
    contigs = {}
    for line in F_contigs.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        contigs[line[1]] = line[0]
    F_contigs.close()
    F_contigs = open(pass_contigs, 'r')
    for line in F_contigs.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        contigs[line[1]] = line[0]
    F_contigs.close()
    return contigs

def cluster_info(cluster_path):
    files = os.listdir(cluster_path)
    contig_id = {}
    for file in files:
        if file[-2:] == "fa":
            for record in SeqIO.parse(cluster_path + file, "fasta"):
                contig_id[str(record.id)] = file.split(".")[0]
    return contig_id


def contig_info(contig_path):
    contigs_seq = {}
    for record in SeqIO.parse(contig_path, "fasta"):
        contigs_seq[str(record.id)] = str(record.seq)
    return contigs_seq


def add_contigs(cluster_path, contigs, contig_id, contigs_seq):
    for key in contigs:
        f_write = open(cluster_path + contig_id[contigs[key]] + ".fa", 'a')
        f_write.write(">" + key + '\n')
        f_write.write(contigs_seq[key] + '\n')
        f_write.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add other types of contigs into the current cluster")
    parser.add_argument("--ensure_contigs", help="Path of Ensure_contigs.txt", required=True, default=None)
    parser.add_argument("--pass_contigs", help="Path of pass_contigs.txt", required=True, default=None)
    parser.add_argument("--cluster_folder",  help="Folder of remaining_cluster.fa", required=True, default=None)
    parser.add_argument("--contig_path", help="Path of othertype_contig.fa", required=True,default=None)
    FLAGS = parser.parse_args()
    contigs = pass_contigs(Ensure_contigs = FLAGS.ensure_contigs, pass_contigs = FLAGS.pass_contigs)
    contig_id = cluster_info(cluster_path = FLAGS.cluster_folder)
    contigs_seq = contig_info(contig_path = FLAGS.contig_folder)
    add_contigs(cluster_path = FLAGS.cluster_path, contigs = contigs, contig_id = contig_id, contigs_seq = contigs_seq)


