import os
from Bio import SeqIO
import argparse
from collections import Counter
import collections


def overlap_contigs(placed_align_path, contigs_fai):
    f = open(contigs_fai, 'r')
    contig_len = {}
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        contig_len[line[0]] = int(line[1])
    f.close()
    contig_cluster = {}
    contig_type = {}
    delete_contigs = set()
    f = open(placed_align_path, 'r')
    for line in f.readlines():
        cluster = set()
        line = line.split('\n')[0]
        line = line.split('\t')
        length = 0
        label = 'all'
        for contig in line:
            if len(contig) > 2:

                ID = contig.split('.')[0]
                if contig_len[ID] >= length:
                    if label != 'b':
                        length = contig_len[ID]
                        rep = ID
                        type = contig.split('.')[1]
                        if type == 'b':
                            label = 'b'
                    else:
                        if contig.split('.')[1] == 'b':
                            length = contig_len[ID]
                            rep = ID
                cluster.add(ID)
                delete_contigs.add(ID)
        contig_cluster[rep] = cluster
        contig_type[rep] = type
    f.close()
    #    print len(contig_cluster)
    return contig_cluster, contig_type, delete_contigs


def obtain_contigs(BEP_cluster_folder, BEP_rep_folder, LEP_cluster_folder, LEP_rep_folder, REP_cluster_folder,
                   REP_rep_folder):
    left_rep = collections.defaultdict(list)
    cluster = collections.defaultdict(list)
    right_rep = collections.defaultdict(list)
    both_rep = collections.defaultdict(list)

    # cluster
    files = os.listdir(LEP_cluster_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(LEP_cluster_folder + file, "fasta"):
            cluster[name].append([record.id, str(record.seq)])

    files = os.listdir(REP_cluster_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(REP_cluster_folder + file, "fasta"):
            cluster[name].append([record.id, str(record.seq)])

    files = os.listdir(BEP_cluster_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(BEP_cluster_folder + file, "fasta"):
            cluster[name].append([record.id, str(record.seq)])

    # rep
    files = os.listdir(LEP_rep_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(LEP_rep_folder + file, "fasta"):
            left_rep[name].append(str(record.seq))

    files = os.listdir(REP_rep_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(REP_rep_folder + file, "fasta"):
            right_rep[name].append(str(record.seq))

    files = os.listdir(BEP_rep_folder)
    for file in files:
        rep = file.split(".")
        if len(rep) == 2:
            name = rep[0]
        elif len(rep) == 3:
            name = rep[0] + "." + rep[1]
        for record in SeqIO.parse(BEP_rep_folder + file, "fasta"):
            both_rep[name].append(str(record.seq))
    return left_rep, right_rep, both_rep, cluster


def generate_placed_contigs(contig_cluster, contig_type, delete_contigs, left_rep, right_rep, both_rep, cluster,
                            LEP_cluster_update_folder, LEP_rep_update_folder, REP_cluster_update_folder,
                            REP_rep_update_folder, BEP_cluster_update_folder, BEP_rep_update_folder):
    for key in right_rep:
        if key not in delete_contigs:
            right_rep_path = open(REP_rep_update_folder + key + '.fa', 'w')
            right_cluster_path = open(REP_cluster_update_folder + key + '.fa', 'w')
            right_rep_path.write('>' + key + '\n')
            right_rep_path.write(right_rep[key][0] + '\n')
            for subcontig in (cluster[key]):
                right_cluster_path.write('>' + subcontig[0] + '\n')
                right_cluster_path.write(subcontig[1] + '\n')
            right_rep_path.close()
            right_cluster_path.close()

    for key in left_rep:
        if key not in delete_contigs:
            left_rep_path = open(LEP_rep_update_folder + key + '.fa', 'w')
            left_cluster_path = open(LEP_cluster_update_folder + key + '.fa', 'w')
            left_rep_path.write('>' + key + '\n')
            left_rep_path.write(left_rep[key][0] + '\n')

            for subcontig in cluster[key]:
                left_cluster_path.write('>' + subcontig[0] + '\n')
                left_cluster_path.write(subcontig[1] + '\n')
            left_rep_path.close()
            left_cluster_path.close()

    for key in both_rep:
        if key not in delete_contigs:
            both_rep_path = open(BEP_rep_update_folder + key + '.fa', 'w')
            both_cluster_path = open(BEP_cluster_update_folder + key + '.fa', 'w')
            both_rep_path.write('>' + key + '\n')
            both_rep_path.write(both_rep[key][0] + '\n')
            for subcontig in cluster[key]:
                both_cluster_path.write('>' + subcontig[0] + '\n')
                both_cluster_path.write(subcontig[1] + '\n')
            both_rep_path.close()
            both_cluster_path.close()

    for contig in contig_cluster:
        #        print contig_type[contig]
        if contig_type[contig] == 'b':
            seq = both_rep[contig][0]
            rep_path = open(BEP_rep_update_folder + contig + '.fa', 'w')
            cluster_path = open(BEP_cluster_update_folder + contig + '.fa', 'w')
        elif contig_type[contig] == 'l':
            seq = left_rep[contig][0]
            rep_path = open(LEP_rep_update_folder + contig + '.fa', 'w')
            cluster_path = open(LEP_cluster_update_folder + contig + '.fa', 'w')
        elif contig_type[contig] == 'r':
            seq = right_rep[contig][0]
            rep_path = open(REP_rep_update_folder + contig + '.fa', 'w')
            cluster_path = open(REP_cluster_update_folder + contig + '.fa', 'w')
        rep_path.write('>' + contig + '\n')
        rep_path.write(seq + '\n')
        for contig_name in (contig_cluster[contig]):
            for subcontig in cluster[contig_name]:
                cluster_path.write('>' + subcontig[0] + '\n')
                cluster_path.write(subcontig[1] + '\n')
        cluster_path.close()
        rep_path.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get the final placed reps and clusters")
    parser.add_argument("--placed_align_path", help="path of placed_aligned.update.tsv", required=True, default=None)
    parser.add_argument("--contigs_fai", help="folder of contig.fa.fai", required=True, default=None)
    parser.add_argument("--BEP_cluster_folder", help="folder of BEP_contigs.fa", required=True, default=None)
    parser.add_argument("--LEP_cluster_folder", help="folder of LEP_contigs.fa", required=True, default=None)
    parser.add_argument("--REP_cluster_folder", help="folder of REP_cluster.fa", required=True, default=None)
    parser.add_argument("--BEP_rep_folder", help="folder of BEP_rep.fa", required=True, default=None)
    parser.add_argument("--LEP_rep_folder", help="folder of LEP_rep.fa", required=True, default=None)
    parser.add_argument("--REP_rep_folder", help="folder of REP_rep.fa", required=True, default=None)
    parser.add_argument("--BEP_cluster_update_folder", help="Save folder of updated BEP_contigs.fa", required=True,
                        default=None)
    parser.add_argument("--LEP_cluster_update_folder", help="Save folder of updated LEP_contigs.fa", required=True,
                        default=None)
    parser.add_argument("--REP_cluster_update_folder", help="Save folder of updated REP_contigs.fa", required=True,
                        default=None)
    parser.add_argument("--BEP_rep_update_folder", help="Save folder of updated BEP_rep.fa", required=True,
                        default=None)
    parser.add_argument("--LEP_rep_update_folder", help="Save folder of updated LEP_rep.fa", required=True,
                        default=None)
    parser.add_argument("--REP_rep_update_folder", help="Save folder of updated REP_rep.fa", required=True,
                        default=None)

    FLAGS = parser.parse_args()

    contig_cluster, contig_type, delete_contigs = overlap_contigs(placed_align_path=FLAGS.placed_align_path,
                                                                  contigs_fai=FLAGS.contigs_fai)
    left_rep, right_rep, both_rep, cluster = obtain_contigs(BEP_cluster_folder=FLAGS.BEP_cluster_folder,
                                                            BEP_rep_folder=FLAGS.BEP_rep_folder,
                                                            LEP_cluster_folder=FLAGS.LEP_cluster_folder,
                                                            LEP_rep_folder=FLAGS.LEP_cluster_folder,
                                                            REP_cluster_folder=FLAGS.REP_cluster_folder,
                                                            REP_rep_folder=FLAGS.REP_rep_folder)
    generate_placed_contigs(contig_cluster, contig_type, delete_contigs, left_rep, right_rep, both_rep, cluster,
                            LEP_cluster_update_folder=FLAGS.LEP_cluster_update_folder,
                            LEP_rep_update_folder=FLAGS.LEP_rep_update_folder,
                            REP_cluster_update_folder=FLAGS.REP_cluster_update_folder,
                            REP_rep_update_folder=FLAGS.REP_rep_update_folder,
                            BEP_cluster_update_folder=FLAGS.BEP_cluster_update_folder,
                            BEP_rep_update_folder=FLAGS.BEP_rep_update_folder)


