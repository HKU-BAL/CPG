import argparse
from Bio import SeqIO



def rep_cluster(seq_path, path_merge_bed, path_rep, path_contig):
    seq = {}
    for record in SeqIO.parse(seq_path, 'fasta'):
        seq[str(record.id)] = [str(record.seq), str(record.seq.reverse_complement())]
    f = open(path_merge_bed, 'r')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        length = int(line[2]) - int(line[1])
        contigs = line[3]
        contigs = contigs.split(',')
        contig_name1 = []
        for contig_info in contigs:
            name = contig_info.split(".")[0].split('_')
            contig_name = name[0] + '_' + name[1]
            sample_name = ''
            for num in range(2, len(name) - 1):
                sample_name += name[num]
                sample_name += '_'
            sample_name += name[len(name) - 1]
            reverse1 = contig_info.split(".")[1].split('_')[1]
            length1 = contig_info.split(".")[1].split('_')[0]
            contig_name1.append([int(length1), contig_info.split(".")[0], contig_name, sample_name, reverse1])
        contig_name = sorted(contig_name1, key=(lambda x: x[0]), reverse=True)
        R_name = contig_name[0][1]
        reverse = contig_name[0][4]
        f_write_R = open(path_rep + R_name + ".fa", 'w')
        if reverse == "+":

            f_write_R.write('>' + R_name + '\n')
            f_write_R.write(seq[R_name][0] + '\n')


        elif reverse == "-":
            f_write_R.write('>' + R_name + '\n')
            f_write_R.write(seq[R_name][1] + '\n')
        f_write_R.close()

        f_write_contig = open(path_contig + R_name + '.fa', 'w')
        for contig_num in range(len(contig_name)):
            if contig_name[contig_num][4] == "+":
                f_write_contig.write('>' + contig_name[contig_num][1] + '\n')
                f_write_contig.write(seq[contig_name[contig_num][1]][0] + '\n')
            elif   contig_name[contig_num][4] == "-":
                f_write_contig.write('>' + contig_name[contig_num][1]+ '\n')
                f_write_contig.write(seq[contig_name[contig_num][1]][1] + '\n')
        f_write_contig.close()
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="get the representatives and sub contigs in each cluster")
    parser.add_argument('--path_merged_bed', help='merged bed path',required=True, default=None)
    parser.add_argument('--path_rep', help='save path of representative ', required=True, default=None)
    parser.add_argument('--path_contig', help='save path of cluster ', required=True, default=None)
    parser.add_argument('--seq_path', help='save path of cluster ', required=True, default=None)
    FLAGS = parser.parse_args()
    rep_cluster(seq_path = FLAGS.seq_path, path_merge_bed = FLAGS.path_merged_bed, path_rep = FLAGS.path_rep, path_contig = FLAGS.path_contig)
