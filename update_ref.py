import os
import collections
from Bio import SeqIO
import argparse

def insertion_points(LEP_folder, REP_folder, contigs_fai):
    f = open(contigs_fai, 'r')
    contig_len = {}
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        contig_len[line[0]] = int(line[1])
    f.close()

    contig_insertion = {}
    files = os.listdir(LEP_folder)
    for file in files:
        f = open(LEP_folder + file, 'r')
        sample_ID = file.split('.txt')[0]
        for line in f.readlines():
            line = line.split('\n')[0]
            line = line.split('\t')
            contig_ID = line[0] + '_' + sample_ID

            length = contig_len[contig_ID]
            if line[3] == "reverse":
                map_len = 101 - int(line[8])
            else:
                map_len = int(line[9]) - 1
            contig_insertion[contig_ID] = [int(line[7]) - 1, 1.0 * (int(line[9]) - int(line[8])) / length, length,map_len]
        f.close()

    files = os.listdir(REP_folder)
    for file in files:
        f = open(REP_folder + file, 'r')
        sample_ID = file.split('.txt')[0]
        for line in f.readlines():
            line = line.split('\n')[0]
            line = line.split('\t')
            contig_ID = line[0] + '_' + sample_ID
            length = contig_len[contig_ID]
            if line[3] == "no_reverse":
                map_len = 101 - int(line[8])
            else:
                map_len = int(line[9]) - 1
            contig_insertion[contig_ID] = [int(line[6]) - 1, 1.0 * (int(line[9]) - int(line[8])) / length, length,map_len]
        f.close()

    return contig_insertion


def obtain_seq(LEP_cluster_folder, REP_cluster_folder, LEP_rep_folder, REP_rep_folder):
    LEP_rep = collections.defaultdict(list)
    LEP_cluster = collections.defaultdict(list)
    REP_rep = collections.defaultdict(list)
    REP_cluster = collections.defaultdict(list)

    # cluster
    files = os.listdir(LEP_cluster_folder)
    for file in files:
        name = file.split(".")[0]
        for record in SeqIO.parse(LEP_cluster_folder + file, "fasta"):
            LEP_cluster[name].append([record.id, str(record.seq)])

    files = os.listdir(REP_cluster_folder)
    for file in files:
        name = file.split(".")[0]
        for record in SeqIO.parse(REP_cluster_folder + file, "fasta"):
            REP_cluster[name].append([record.id, str(record.seq)])

    # rep
    files = os.listdir(LEP_rep_folder)
    for file in files:
        name = file.split(".")[0]
        for record in SeqIO.parse(LEP_rep_folder + file, "fasta"):
            LEP_rep[name].append(str(record.seq))

    files = os.listdir(REP_rep_folder)
    for file in files:
        name = file.split(".")[0]
        for record in SeqIO.parse(REP_rep_folder + file, "fasta"):
            REP_rep[name].append(str(record.seq))
    return LEP_cluster, REP_cluster, LEP_rep, REP_rep


def move_SEP_contigs(file, LEP_cluster, REP_cluster, LEP_add, REP_add, delete_LEP, delete_REP):
    f = open(file, 'r')

    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')

        rep = line[len(line) - 1].split('.')
        if rep[len(rep) - 1] == 'l':
            for number1 in range(len(LEP_cluster[line[0]])):
                LEP_add[rep[0]].append([LEP_cluster[line[0]][number1][0], LEP_cluster[line[0]][number1][1]])

            for number1 in range(len(REP_cluster[line[1]])):
                LEP_add[rep[0]].append([REP_cluster[line[1]][number1][0], REP_cluster[line[1]][number1][1]])

        elif rep[len(rep) - 1] == 'r':
            for number1 in range(len(LEP_cluster[line[0]])):
                REP_add[rep[0]].append([LEP_cluster[line[0]][number1][0], LEP_cluster[line[0]][number1][1]])

            for number1 in range(len(REP_cluster[line[1]])):
                REP_add[rep[0]].append([REP_cluster[line[1]][number1][0], REP_cluster[line[1]][number1][1]])

        delete_LEP.add(line[0])
        delete_REP.add(line[1])
    f.close()
    return LEP_add, REP_add, delete_LEP, delete_REP


def overlap_SEP(file, contig_insert_pos):
    f = open(file, 'r')
    move_LEP = set()
    move_REP = set()
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        seq1 = line[len(line) - 2]  # LEP
        seq2 = line[len(line) - 1]  # REP
        contig_ID = seq1 + "." + seq2
        if int(contig_insert_pos[seq1][0]) > int(contig_insert_pos[seq2][0]):
            if contig_insert_pos[seq1][1] < contig_insert_pos[seq2][1]:
                move_LEP.add(contig_ID)
            elif contig_insert_pos[seq1][1] > contig_insert_pos[seq2][1]:
                move_REP.add(contig_ID)
            else:
                if contig_insert_pos[seq1][2] < contig_insert_pos[seq2][2]:
                    move_LEP.add(contig_ID)
                elif contig_insert_pos[seq1][2] > contig_insert_pos[seq2][2]:
                    move_REP.add(contig_ID)
                else:
                    if contig_insert_pos[seq1][3] < contig_insert_pos[seq2][3]:
                        move_LEP.add(contig_ID)
                    elif contig_insert_pos[seq1][3] > contig_insert_pos[seq2][3]:
                        move_REP.add(contig_ID)
    f.close()
    return move_LEP, move_REP


def check_BEP_contigs(file, move_LEP, move_REP, LEP_cluster, REP_cluster, delete_LEP, delete_REP, LEP_add, REP_add,LEP_rep, REP_rep, merged_contig_folder):
    f = open(file, 'r')
    BEP_add = collections.defaultdict(list)
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        LEP_ID = line[len(line) - 2]
        REP_ID = line[len(line) - 1]
        rep = LEP_ID + "." + REP_ID
        for record in SeqIO.parse(merged_contig_folder + rep + "/" + rep + '.fa', "fasta"):
            sequence = str(record.seq)
        if rep in move_LEP:
            LEP_rep[rep] = [sequence]
            for number1 in range(len(LEP_cluster[line[0]])):
                LEP_add[rep].append([LEP_cluster[line[0]][number1][0], LEP_cluster[line[0]][number1][1]])
            for number1 in range(len(REP_cluster[line[1]])):
                LEP_add[rep].append([REP_cluster[line[1]][number1][0], REP_cluster[line[1]][number1][1]])
        elif rep in move_REP:
            REP_rep[rep] = [sequence]
            for number1 in range(len(LEP_cluster[line[0]])):
                REP_add[rep].append([LEP_cluster[line[0]][number1][0], LEP_cluster[line[0]][number1][1]])
            for number1 in range(len(REP_cluster[line[1]])):
                REP_add[rep].append([REP_cluster[line[1]][number1][0], REP_cluster[line[1]][number1][1]])
        else:
            for number1 in range(len(LEP_cluster[line[0]])):
                BEP_add[rep].append([LEP_cluster[line[0]][number1][0], LEP_cluster[line[0]][number1][1]])
            for number1 in range(len(REP_cluster[line[1]])):
                BEP_add[rep].append([REP_cluster[line[1]][number1][0], REP_cluster[line[1]][number1][1]])

        delete_LEP.add(line[0])
        delete_REP.add(line[1])
    f.close()
    return BEP_add, LEP_rep, REP_rep, LEP_add, REP_add, delete_LEP, delete_REP


def reobtain_rep_cluster(LEP_cluster, REP_cluster, LEP_rep, REP_rep, BEP_add, LEP_add, REP_add, delete_LEP, delete_REP,
                         REP_rep_update_folder, REP_cluster_update_folder, LEP_rep_update_folder,
                         LEP_cluster_update_folder, BEP_rep_folder, BEP_cluster_folder, merged_contig_folder):
    for key in REP_rep:
        if key not in delete_REP:
            right_rep_path = open(REP_rep_update_folder + '/' + key + '.fa', 'w')
            right_cluster_path = open(REP_cluster_update_folder + '/' + key + '.fa', 'w')
            right_rep_path.write('>' + key + '\n')
            right_rep_path.write(REP_rep[key][0] + '\n')

            for number in range(len(REP_cluster[key])):
                right_cluster_path.write('>' + REP_cluster[key][number][0] + '\n')
                right_cluster_path.write(REP_cluster[key][number][1] + '\n')
            right_rep_path.close()
            right_cluster_path.close()

    for key in LEP_rep:
        if key not in delete_LEP:
            left_rep_path = open(LEP_rep_update_folder + '/' + key + '.fa', 'w')
            left_cluster_path = open(LEP_cluster_update_folder + '/' + key + '.fa', 'w')
            left_rep_path.write('>' + key + '\n')
            left_rep_path.write(LEP_rep[key][0] + '\n')
            for number in range(len(LEP_cluster[key])):
                left_cluster_path.write('>' + LEP_cluster[key][number][0] + '\n')
                left_cluster_path.write(LEP_cluster[key][number][1] + '\n')
            left_cluster_path.close()
            left_rep_path.close()

    for key in REP_add:
        right_rep_path = open(REP_rep_update_folder + '/' + key + '.fa', 'w')
        right_cluster_path = open(REP_cluster_update_folder + '/' + key + '.fa', 'w')
        right_rep_path.write('>' + key + '\n')
        right_rep_path.write(REP_rep[key][0] + '\n')
        for number in range(len(REP_add[key])):
            right_cluster_path.write('>' + REP_add[key][number][0] + '\n')
            right_cluster_path.write(REP_add[key][number][1] + '\n')
        right_rep_path.close()
        right_cluster_path.close()

    for key in LEP_add:
        left_rep_path = open(LEP_rep_update_folder + "/" + key + '.fa', 'w')
        left_cluster_path = open(LEP_cluster_update_folder + "/" + key + '.fa', 'w')
        left_rep_path.write('>' + key + '\n')
        left_rep_path.write(LEP_rep[key][0] + '\n')
        for number in range(len(LEP_add[key])):
            left_cluster_path.write('>' + LEP_add[key][number][0] + '\n')
            left_cluster_path.write(LEP_add[key][number][1] + '\n')
        left_rep_path.close()
        left_cluster_path.close()

    # add 2EP contigs
    for key in BEP_add:
        both_rep_path = open(BEP_rep_folder + '/' + key + '.fa', 'w')
        both_cluster_path = open(BEP_cluster_folder + '/' + key + '.fa', 'w')

        for record in SeqIO.parse(merged_contig_folder + key + "/" + key + '.fa', "fasta"):
            both_rep_path.write('>' + key + '\n')
            both_rep_path.write(str(record.seq) + '\n')

        for number in range(len(BEP_add[key])):
            both_cluster_path.write('>' + BEP_add[key][number][0] + '\n')
            both_cluster_path.write(BEP_add[key][number][1] + '\n')
        both_cluster_path.close()
        both_rep_path.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Update the placed representatives and corresponding clusters.")
    parser.add_argument("--LEP_folder",help="the folder of LEP contigs' alignment information (step2.files from point6)", required=True, default=None)
    parser.add_argument("--REP_folder", help="fthe folder of REP contigs' alignment information (step2.files from point6)", required=True, default=None)
    parser.add_argument("--contigs_fai", help="folder of contig.fa.fai", required=True, default=None)
    parser.add_argument("--LEP_cluster_folder", help="folder of LEP_contigs.fa", required=True, default=None)
    parser.add_argument("--REP_cluster_folder", help="folder of REP_cluster.fa", required=True, default=None)
    parser.add_argument("--LEP_rep_folder", help="folder of LEP_rep.fa", required=True, default=None)
    parser.add_argument("--REP_rep_folder", help="folder of REP_rep.fa", required=True, default=None)
    parser.add_argument("--align_folder", help="folder of alignment result", required=True, default=None)
    parser.add_argument("--merged_contig_folder", help="Path of merge.fa generated by PopINS", required=True,default=None)
    parser.add_argument("--LEP_rep_update_folder", help="Save folder of updated LEP_rep.fa", required=True, default=None)
    parser.add_argument("--LEP_cluster_update_folder", help="Save folder of updated LEP_contigs.fa", required=True,default=None)
    parser.add_argument("--REP_rep_update_folder", help="Save folder of updated REP_rep.fa", required=True,default=None)
    parser.add_argument("--REP_cluster_update_folder", help="Save folder of updated REP_contigs.fa", required=True, default=None)
    parser.add_argument("--BEP_rep_folder", help="folder of REP_rep.fa", required=True, default=None)
    parser.add_argument("--BEP_cluster_folder", help="folder of BEP_contigs.fa", required=True, default=None)
    FLAGS = parser.parse_args()
    
    LEP_cluster, REP_cluster, LEP_rep, REP_rep = obtain_seq(LEP_cluster_folder=FLAGS.LEP_cluster_folder,
                                                            REP_cluster_folder=FLAGS.REP_cluster_folder,
                                                            LEP_rep_folder=FLAGS.LEP_rep_folder,
                                                            REP_rep_folder=FLAGS.REP_rep_folder)
    alignment_files = ["final_identity.txt", "final_contained.txt", "final_part.txt"]
    LEP_add = collections.defaultdict(list)
    REP_add = collections.defaultdict(list)
    delete_LEP = set()
    delete_REP = set()
    for file in alignment_files:
        file_path = FLAGS.align_folder + file
        LEP_add, REP_add, delete_LEP, delete_REP = move_SEP_contigs(file_path, LEP_cluster, REP_cluster, LEP_add,
                                                                    REP_add, delete_LEP, delete_REP)
    file = FLAGS.align_folder + "final_overlap.txt"
    contig_insertion = insertion_points(LEP_folder=FLAGS.LEP_folder, REP_folder=FLAGS.REP_folder,
                                        contigs_fai=FLAGS.contigs_fai)
    move_LEP, move_REP = overlap_SEP(file, contig_insertion)
    BEP_add, LEP_rep, REP_rep, LEP_add, REP_add, delete_LEP, delete_REP = check_BEP_contigs(file, move_LEP, move_REP,
                                                                                            LEP_cluster, REP_cluster,
                                                                                            delete_LEP, delete_REP,
                                                                                            LEP_add, REP_add, LEP_rep,
                                                                                            REP_rep,
                                                                                            merged_contig_folder=FLAGS.merged_contig_folder)
    reobtain_rep_cluster(LEP_cluster, REP_cluster, LEP_rep, REP_rep, BEP_add, LEP_add, REP_add, delete_LEP, delete_REP,
                         REP_rep_update_folder=FLAGS.REP_rep_update_folder,
                         REP_cluster_update_folder=FLAGS.REP_cluster_update_folder,
                         LEP_rep_update_folder=FLAGS.LEP_rep_update_folder,
                         LEP_cluster_update_folder=FLAGS.LEP_cluster_update_folder, BEP_rep_folder=FLAGS.BEP_rep_folder,
                         BEP_cluster_folder=FLAGS.BEP_cluster_folder, merged_contig_folder=FLAGS.merged_contig_folder)
