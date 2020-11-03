import os
import argparse

def alignment_priority(LEP_bed, REP_bed, Identity_path, Contained_path, Overlap_path, Part_align_path):
    f = open(LEP_bed, 'r')
    reverse_info = {}
    pos_begin = {}
    pos_end = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        seq_name = line[3].split(".")[0]
        reverse = line[4]
        reverse_info[seq_name] = reverse
        pos_begin[seq_name] = line[1]
        pos_end[seq_name] = line[2]
    f.close()

    f = open(REP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        seq_name = line[3].split(".")[0]
        reverse = line[4]
        reverse_info[seq_name] = reverse
        pos_begin[seq_name] = line[1]
        pos_end[seq_name] = line[2]
    f.close()

    f = open(Identity_path, 'r')
    left_identity = set()
    right_identity = set()
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        left_identity.add(left_contig)
        right_identity.add(right_contig)
    f.close()

    f_contain = open("/dev/shm/contained1.txt", 'w')
    f = open(Contained_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        if left_contig not in left_identity and right_contig not in right_identity:
            for num in range(len(line) - 1):
                f_contain.write(line[num] + '\t')
            f_contain.write(line[len(line) - 1] + '\n')
        left_identity.add(left_contig)
        right_identity.add(right_contig)
    f_contain.close()
    f.close()

    f_overlap = open("/dev/shm/overlap1.txt", 'w')
    f = open(Overlap_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        if reverse_info[left_contig] == reverse_info[right_contig] and pos_begin[left_contig] <= pos_end[right_contig]:
            if left_contig not in left_identity and right_contig not in right_identity:
                for num in range(len(line) - 1):
                    f_overlap.write(line[num] + '\t')
                f_overlap.write(line[len(line) - 1] + '\n')
            left_identity.add(left_contig)
            right_identity.add(right_contig)
    f_overlap.close()
    f.close()

    f_part = open("/dev/shm/part1.txt", 'w')
    f = open(Part_align_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        if left_contig not in left_identity and right_contig not in right_identity:
            for num in range(len(line) - 1):
                f_part.write(line[num] + '\t')
            f_part.write(line[len(line) - 1] + '\n')
    f.close()
    f_part.close()

def redefine_alignment_result(Identity_path, save_folder):
    f = open(Identity_path, 'r')

    f_identity = open(save_folder+"final_identity.txt", 'w')
    right = {}
    left = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]
        left_length = float(line[7])
        right_length = float(line[8])

        if left_contig not in left:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                if right[right_contig][2] < identity:
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,
                                                    right_length]
                elif right[right_contig][2] == identity and right[right_contig][3] < cov1 and right[right_contig][4] < cov2:
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,
                                                    right_length]
        else:
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2:
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
    f.close()

    for key in left:
        f_identity.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_identity.write(str(left[key][0]) + '.l' + '\n')
        else:
            f_identity.write(str(left[key][1]) + '.r' + '\n')
    f_identity.close()

    f_contain = open(save_folder+"final_contained.txt", 'w')
    right = {}
    left = {}
    f = open("/dev/shm/contained1.txt", 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]
        left_length = float(line[7])
        right_length = float(line[8])

        if left_contig not in left:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                if right[right_contig][2] < identity:
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,
                                                    right_length]
                elif right[right_contig][2] == identity and (
                        right[right_contig][3] < cov1 and right[right_contig][4] < cov2):
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]

        else:
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and (left[left_contig][3] < cov1 and left[left_contig][4] and cov2):
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
    f.close()

    for key in left:
        f_contain.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_contain.write(str(left[key][0]) + '.l' + '\n')
        else:
            f_contain.write(str(left[key][1]) + '.r' + '\n')
    f_contain.close()


    f_overlap = open(save_folder +"final_overlap.txt", 'w')
    right = {}
    left = {}
    f = open("/dev/shm/overlap1.txt", 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]

        if left_contig not in left:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2]
            else:
                if right[right_contig][2] < identity:
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2]
                elif right[right_contig][2] == identity and (
                        right[right_contig][3] < cov1 and right[right_contig][4] < cov2):
                    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2]

        else:
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2]
            if left[left_contig][2] < identity:
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2]
            elif left[left_contig][2] == identity and (left[left_contig][3] < cov1 and left[left_contig][4] < cov2):
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2]
    f.close()

    for key in left:
        f_overlap.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\n')
    f_overlap.close()

    f_part = open(save_folder +"f_part.txt", 'w')
    right = {}
    left = {}
    f = open("/dev/shm/part1.txt", 'r')

    length = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]
        left_length = float(line[7])
        right_length = float(line[8])
        length[left_contig] = [left_length, right_length]
        if left_contig not in left:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                if right[right_contig][2] < identity:
                    left[left_length] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                elif right[right_contig][2] == identity and (
                        right[right_contig][3] < cov1 and right[right_contig][4] < cov2):

                    left[left_length] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        else:
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and (left[left_contig][3] < cov1 and left[left_contig][4] < cov2):

                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]

    for key in left:
        f_part.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\n')
    f_part.close()
    f.close()

    if os.path.exists("/dev/shm/overlap1.txt"):
        os.remove("/dev/shm/overlap1.txt")
    if os.path.exists("/dev/shm/contained1.txt"):
        os.remove("/dev/shm/contained1.txt")
    if os.path.exists("/dev/shm/part1.txt"):
        os.remove("/dev/shm/part1.txt")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Based on the priority of alignment, re-obtain the four types of alignment result.")
    parser.add_argument("--LEP_bed", help="Path of LEP_contigs.bed", required=True, default=None)
    parser.add_argument("--REP_bed", help="Path of REP_contigs.bed", required=True, default=None)
    parser.add_argument("--Identity_path", help="Path of identity.coords", required=True, default=None)
    parser.add_argument("--Contained_path", help="Path of contained.coords", required=True, default=None)
    parser.add_argument("--Overlap_path", help="Path of overlap.coords", required=True, default=None)
    parser.add_argument("--Part_align_path", help="Path of part.coords", required=True, default=None)
    parser.add_argument("--save_folder", help="Save the four new alignment results", required=True, default=None)
    FLAGS = parser.parse_args()
    alignment_priority(LEP_bed = FLAGS.LEP_bed, REP_bed= FLAGS.REP_bed, Identity_path= FLAGS.Identity_path, Contained_path = FLAGS.Contained_path, Overlap_path = FLAGS.Overlap_path, Part_align_path = FLAGS.Part_align_path)
    redefine_alignment_result(Identity_path = FLAGS.Identity_path, save_folder = FLAGS.save_folder)
