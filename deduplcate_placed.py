import os
import collections
import argparse
import networkx as nx


def best_hit_alignment(alignment_path, identity, coverage):
    f_write = open("/dev/shm/placed_filter.tsv", 'w')
    f = open(alignment_path, 'r')
    record = 2
    hit_contigs = set()
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')

        cov1 = 1.0 * (float(line[5]) - int(line[12])) / float(line[4])
        cov2 = 1.0 * (float(line[5]) - int(line[12])) / float(line[3])
        if line[0] not in hit_contigs:
            if line[0] == line[1]:
                record = 1
            else:
                record = 2
                hit_contigs.add(line[0])
                if float(line[2]) >= identity and (cov1 >= coverage or cov2 >= coverage):
                    for num in range(len(line) - 1):
                        f_write.write(line[num] + '\t')
                    f_write.write(line[len(line) - 1] + '\n')
        elif line[0] in hit_contigs and record == 1:
            record = 2
            hit_contigs.add(line[0])
            if float(line[2]) >= identity and (cov1 >= coverage or cov2 >= coverage):
                for num in range(len(line) - 1):
                    f_write.write(line[num] + '\t')
                f_write.write(line[len(line) - 1] + '\n')
    f_write.close()
    f.close()


def placement_pos(BEP_bed, LEP_bed, REP_bed):
    f = open(BEP_bed, 'r')
    contig_pos = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split(".")[0]
        start = int(line[1])
        end = int(line[2])
        contig_pos[name] = [str(line[0]), start, end, '.b']
    #        pos.add(name)
    f.close()

    f = open(LEP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split(".")[0]
        start = int(line[1])
        end = int(line[2])
        contig_pos[name] = [str(line[0]), start, end, '.l']
    f.close()

    f = open(REP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split(".")[0]
        start = int(line[1])
        end = int(line[2])
        contig_pos[name] = [str(line[0]), start, end, '.r']
    f.close()

    return contig_pos


def judge_distance(distance, contig_pos, pass_alignment):
    obtain_contigs = collections.defaultdict(list)
    map_contigs = []

    f = open('/mnt/bal19/qhli/CPG1/step3/placed/final_align/placed_highqual.tsv', 'r')
    f_write = open(pass_alignment, 'w')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        q1 = line[0]
        q2 = line[1]
        if q1 not in contig_pos:
            name1 = q1.split('.')[0]
            name2 = q1.split('.')[1]
            query1_chr = contig_pos[name1][0]
            query1_begin = contig_pos[name1][1]
            query1_end = contig_pos[name2][1]
            orentation1 = '.b'

        else:
            query1_chr = contig_pos[q1][0]
            query1_begin = contig_pos[q1][1]
            query1_end = contig_pos[q1][2]
            orentation1 = contig_pos[q1][3]

        if q2 not in contig_pos:
            name1 = q2.split('.')[0]
            name2 = q2.split('.')[1]
            query2_chr = contig_pos[name1][0]
            query2_begin = contig_pos[name1][1]
            query2_end = contig_pos[name2][1]
            orentation2 = '.b'
        else:
            query2_chr = contig_pos[q2][0]
            query2_begin = contig_pos[q2][1]
            query2_end = contig_pos[q2][2]
            orentation2 = contig_pos[q2][3]

        if query2_chr == query1_chr:
            if query2_begin - query1_end < distance and query1_begin - query2_end < distance:
                if q2 not in obtain_contigs:
                    map_contigs.append([q1 + orentation1, q2 + orentation2])
                    obtain_contigs[q1].append(q2)
                else:
                    if q1 not in obtain_contigs[q2]:
                        map_contigs.append([q1 + orentation2, q2 + orentation2])

                        obtain_contigs[q1].append(q2)

    f.close()
    G = nx.Graph()
    G.add_nodes_from(sum(map_contigs, []))
    info = [[(s[i], s[i + 1]) for i in range(len(s) - 1)] for s in map_contigs]
    for i in info:
        G.add_edges_from(i)
    for i in nx.connected_components(G):
        information = ''
        for z in i:
            information += z + '\t'
        f_write.write(information + '\n')  # for i in nx.connected_components(G)]
    os.remove("/dev/shm/placed_filter.tsv")
    f_write.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter out alignment results.")
    parser.add_argument("--alignment_path", help="path of filtered_align.tsv", required=True, default=None)
    parser.add_argument("--BEP_bed", help="path of BEP_contigs.bed ", required=True, default=None)
    parser.add_argument("--LEP_bed", help="folder of LEP_contigs.bed", required=True, default=None)
    parser.add_argument("--REP_bed", help="folder of REP_contigs.bed", required=True, default=None)
    parser.add_argument("--distance", help="the distance between the placement locations of two contigs",
                        required=False, default=2000)
    parser.add_argument("--identity", help="identity cutoff", required=False, default=98)
    parser.add_argument("--coverage", help="coverage cutoff", required=False, default=0.95)
    parser.add_argument("--pass_alignment", help="path of pass alignment.tsv", required=True, default=None)
    FLAGS = parser.parse_args()
    best_hit_alignment(alignment_path=FLAGS.alignment_path, identity=FLAGS.identity, coverage=FLAGS.coverage)
    contig_pos = placement_pos(BEP_bed=FLAGS.BEP_bed, LEP_bed=FLAGS.LEP_bed, REP_bed=FLAGS.REP_bed)
    judge_distance(distance=FLAGS.distance, contig_pos=contig_pos, pass_alignment=FLAGS.pass_alignment)


