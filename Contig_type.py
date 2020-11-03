import os
import argparse

def chr_id_map(ref_name_id):
    f = open(ref_name_id, 'r')
    chr_nc = {}
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        chr_nc[line[0]] = line[1]
    return chr_nc

def single_placed(chr_nc_id, alignment_info, LEP_contigs , REP_contigs):
    ref_set = set()
    for i in range(1, 23):
        ref_set.add(str(i))
    ref_set.add('X')
    ref_set.add('Y')

    for root, dirs, files in os.walk(alignment_info):
        root1 = root.split("/")
        sample_name = root1[len(root1) - 1]
        if sample_name != '':
            f_write_left = open(LEP_contigs + sample_name + '.txt', 'w')
            f_write_right = open(REP_contigs + sample_name + '.txt', 'w')

            for file in files:
                f = open(root + '/' + file, 'r')

                contig_name = file.split('.')[0]
                contig_start = []
                contig_end = []
                ref_start = []
                ref_end = []
                contig_reverse = []
                for line in f.readlines():
                    line = line.split(" ")

                    if len(line) == 4:
                        ref = line[0].split('>')[1]
                        ref = ref.split(':')
                        ref_name = ref[0]
                        if ref_name in ref_set:
                            nc_id = chr_nc_id[ref_name]
                            number = ref[1].split('-')
                            ref_pos_begin = int(number[0])



                    if len(line) == 7:
                        if ref_name not in ref_set:
                            break
                        else:
                            if int(line[2]) < int(line[3]) and int(line[2]) <= 5:

                                contig_start.append(int(line[2]))
                                contig_end.append(int(line[3]))
                                contig_reverse.append(0)

                                ref_start.append(int(line[0]))
                                ref_end.append(int(line[1]))

                            elif int(line[3]) <= 5:

                                contig_start.append(int(line[3]))
                                contig_end.append(int(line[2]))
                                contig_reverse.append(1)
                                ref_start.append(int(line[0]))
                                ref_end.append(int(line[1]))

                if len(contig_start) > 1:
                    if contig_reverse[0] == 0:
                        index_num = ref_start.index(min(ref_start))
                        leftmost_ref_start = min(ref_start)
                        leftmost_contig_start = contig_start[index_num]
                        leftmost_contig_end = contig_end[index_num]
                        length = abs(leftmost_contig_end - leftmost_contig_start)
                        insertion_point = ref_pos_begin + leftmost_ref_start - leftmost_contig_start + 1

                        f_write_left.write(
                            contig_name + '\t' + ref_name + '\t' + str(insertion_point) + '\t' + "no_reverse" + '\t' + str(
                                length) + '\t' + nc_id + '\t' + str(
                                ref_pos_begin + ref_start[index_num]) + '\t' + str(
                                ref_pos_begin + ref_end[index_num]) + '\t' + str(contig_start[index_num]) + '\t' + str(
                                contig_end[index_num]) + '\t' + 'more' + '\n')
                    else:
                        index_num = ref_end.index(max(ref_end))

                        leftmost_ref_end = ref_end[index_num]
                        leftmost_contig_start = contig_start[index_num]
                        leftmost_contig_end = contig_end[index_num]
                        length = abs(leftmost_contig_end - leftmost_contig_start)
                        insertion_point = ref_pos_begin + leftmost_ref_end + leftmost_contig_start - 1

                        f_write_right.write(
                            contig_name + '\t' + ref_name + '\t' + str(insertion_point) + '\t' + "reverse" + '\t' + str(
                                length) + '\t' + nc_id + '\t' + str(
                                ref_pos_begin + ref_start[index_num]) + '\t' + str(
                                ref_pos_begin + ref_end[index_num]) + '\t' + str(contig_start[index_num]) + '\t' + str(
                                contig_end[index_num]) + '\t' + 'more' + '\n')
                elif len(contig_start) == 1:
                    if contig_reverse[0] == 0:
                        insertion_point = ref_pos_begin + ref_start[0] - contig_start[0] + 1
                        length = abs(contig_end[0] - contig_start[0])

                        f_write_left.write(
                            contig_name + '\t' + ref_name + '\t' + str(insertion_point) + '\t' + "no_reverse" + '\t' + str(
                                length) + '\t' + nc_id + '\t' + str(
                                ref_pos_begin + ref_start[0]) + '\t' + str(ref_pos_begin + ref_end[0]) + '\t' + str(
                                contig_start[0]) + '\t' + str(
                                contig_end[0]) + '\n')
                    else:
                        insertion_point = ref_pos_begin + ref_end[0] + contig_start[0] - 1
                        length = abs(contig_end[0] - contig_start[0])

                        f_write_right.write(
                            contig_name + '\t' + ref_name + '\t' + str(insertion_point) + '\t' + "reverse" + '\t' + str(
                                length) + '\t' + nc_id + '\t' + str(
                                ref_pos_begin + ref_start[0]) + '\t' + str(ref_pos_begin + ref_end[0]) + '\t' + str(
                                contig_start[0]) + '\t' + str(
                                contig_end[0]) + '\n')
            f_write_left.close()
            f_write_right.close()


def both_placed(LEP_contigs , REP_contigs, BEP_contigs ,BEP_contigs_all):
    files = os.listdir(LEP_contigs)
    for file in files:
        # print file
        left_contigs_info = {}
        f_write_both = open(BEP_contigs + file, 'w')
        f_write_both_all = open(BEP_contigs_all + file, 'w')
        f = open(LEP_contigs + file)
        for line in f.readlines():
            line = line.split('\t')
            left_contigs_info[line[0]] = ([line[1], line[2], line[3], line[4]])

        right_contigs_info = {}
        f = open(REP_contigs+ file)
        for line in f.readlines():
            line = line.split('\t')
            right_contigs_info[line[0]] = ([line[1], line[2], line[3], line[4]])

        for key in left_contigs_info:
            if key in right_contigs_info:
                left_chr = left_contigs_info[key][0]
                right_chr = right_contigs_info[key][0]
                f_write_both_all.write(key+'\n')
                if left_chr == right_chr:
                    left_reverse = left_contigs_info[key][2]
                    right_reverse = right_contigs_info[key][2]
                    if left_reverse == right_reverse:
                        if left_reverse == 'reverse':
                            if float(right_contigs_info[key][1]) <= float(left_contigs_info[key][1]):
                                left_length = left_contigs_info[key][3]
                                right_length = right_contigs_info[key][3]
                                left_insertion = left_contigs_info[key][1]
                                right_insertion = right_contigs_info[key][1]
                                f_write_both.write(
                                    key + '\t' + left_chr + '\t' + left_insertion + '\t' + left_length + '\t' + right_insertion + '\t' + right_length + '\t' + left_reverse + '\n')

                        elif left_reverse == 'no_reverse':
                            if float(left_contigs_info[key][1]) <= float(right_contigs_info[key][1]):
                                left_length = left_contigs_info[key][3]
                                right_length = right_contigs_info[key][3]
                                left_insertion = left_contigs_info[key][1]
                                right_insertion = right_contigs_info[key][1]
                                f_write_both.write(
                                    key + '\t' + left_chr + '\t' + left_insertion + '\t' + left_length + '\t' + right_insertion + '\t' + right_length + '\t' + left_reverse + '\n')
        f_write_both.close()
        f_write_both_all.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Determine BEP/SEP/Unplaced contigs")
    parser.add_argument("--ref_name_id", help="GRCh38 name to ID ", required=True, default=None)
    parser.add_argument("--alignment_info", help="the path of filtered_info.delta", required=True, default=None)
    parser.add_argument("--LEP_contigs",  help="save path of LEP contigs placed information ", required=True, default=None)
    parser.add_argument("--REP_contigs", help="save path of REP contigs placed information ", required=True,default=None)
    parser.add_argument("--BEP_contigs", help="save path of BEP contigs placed information ", required=True, default=None)
    parser.add_argument("--BEP_contigs_all", help="save path of contigs that left and right ends are placed on reference ", required=True, default=None)
    FLAGS = parser.parse_args()
    chr_nc_id = chr_id_map( ref_name_id = FLAGS.ref_name_id)
    single_placed(chr_nc_id = chr_nc_id, alignment_info = FLAGS.alignment_info, LEP_contigs = FLAGS.LEP_contigs, REP_contigs = FLAGS.REP_contigs)
    both_placed(LEP_contigs = FLAGS.LEP_contigs, REP_contigs = FLAGS.REP_contigs, BEP_contigs = FLAGS.BEP_contigs, BEP_contigs_all = FLAGS.BEP_contigs_all)




