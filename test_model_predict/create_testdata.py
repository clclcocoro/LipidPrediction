#!/usr/bin/env python

import sys
import random

Usage = """
Usage: create_testdata.py <interval> <sequence_length> [--randomscore]

    <interval>        integer that is bigger than 8
    <sequence_length> integer that is bigger than 9
    --randomscore     score is sampled from discrete random distribution
                      binding residue                             randint(3, 5)
                      the both sides residues of binding residue  randint(-10, -8)
                      non-binding residue                         randint(-2, 0)
"""

     
pssm_format = """
Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{}
                      K         Lambda
Standard Ungapped    0.1362     0.3170
Standard Gapped      0.0520     0.2670
PSI Ungapped         0.1704     0.3120
PSI Gapped           0.0520     0.2670
"""


# start is the first binding residue index.
# interval is the number of residues between binding residues.
def generate_pssm(start, sequence_length, interval, random_flag=False):
    interval += 1 # Modify interval.
    pssm = []
    if start == 0:
        neighbor_index = [interval-1, start+1]
    elif start == interval:
        neighbor_index = [start-1, 0]
    else:
        neighbor_index = [start-1, start+1]
    for i in xrange(sequence_length):
        if i % interval == start:
            if random_flag:
                pssm.append([random.randint(3, 5) for i in xrange(20)])
            else:
                pssm.append([2]*20)
        elif i % interval == neighbor_index[0] or i % interval == neighbor_index[1]:
            if random_flag:
                pssm.append([random.randint(-10, -8) for i in xrange(20)])
            else:
                pssm.append([1]*20)
        else:
            if random_flag:
                pssm.append([random.randint(-2, 0) for i in xrange(20)])
            else:
                pssm.append([-1]*20)
    return pssm


def generate_sequence(amino_acid_list, start_aa_index, sequence_length):
    sequence = ""
    for i in xrange(sequence_length):
        j = i / 10
        sequence += amino_acid_list[start_aa_index+j]
    return sequence


def generate_binding_record(sequence_length, start_bindres_index, interval):
    bind_record = ""
    if start_bindres_index != 0:
        bind_record += "0" * start_bindres_index
    for i in xrange(sequence_length-start_bindres_index):
        if i % (interval+1) == 0:
            bind_record += "1"
        else:
            bind_record += "0"
    return bind_record


def write_to_bindres_file(bindres_dir, filename, proteinid, sequence, bind_record):
    with open("{}/{}".format(bindres_dir, filename), 'w') as fp:
        fp.write(">{}\n{}\n{}".format(proteinid, sequence, bind_record))


def generate_formatted_pssm(pssm, sequence):
    formatted_pssm = ""
    for i, row in enumerate(pssm):
        amino_acid = sequence[i]
        str_row = ""
        for ele in row:
            str_row += "{:3d}".format(ele)
        formatted_pssm +=" {:4d} {}  {}    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  0.00 0.00\n".format(i+1, amino_acid, str_row)
    return formatted_pssm


def write_to_pssm_file(pssm_dir, filename, pssm, sequence):
    formatted_pssm = generate_formatted_pssm(pssm, sequence)
    with open("{}/{}".format(pssm_dir, filename), 'w') as fp:
        fp.write(pssm_format.format(formatted_pssm))


# 'H' if binding residue. Both side residues of binding residue are 'E'.
def generate_secondary_structure(bind_record, sequence_length):
    secondary_structure = ['-']*sequence_length
    bind_indexes = []
    for i, bind in enumerate(bind_record):
        if bind == '1':
            bind_indexes.append(i)
    for bind_index in bind_indexes:
        if bind_index != 0:
            secondary_structure[bind_index-1] = 'E'
        if bind_index != sequence_length-1:
            secondary_structure[bind_index+1] = 'E'
    for bind_index in bind_indexes:
        secondary_structure[bind_index] = 'H'
    return ''.join(secondary_structure)




def write_to_secondary_structure_file(secondary_structure_dir, filename, proteinid, secondary_structure, sequence):
    with open("{}/{}".format(secondary_structure_dir, filename), 'w') as fp:
        fp.write(">{}\n{}\n{}".format(proteinid, sequence, secondary_structure))


def generate_proteins_data(number_of_proteins, bindres_dir, pssm_dir, secondary_structure_dir, sequence_length, interval, random_flag, positive=True):
    data_paths = ""
    for i in xrange(number_of_proteins):
        if positive:
            proteinid = "PPPP:{}".format(i)
        else:
            proteinid = "NNNN:{}".format(i)
        bindres_filename = "{}.bindres".format(proteinid)
        pssm_filename = "{}.pssm".format(proteinid)
        secondary_structure_filename = "{}.secondary_structure.txt".format(proteinid)

        sequence = generate_sequence(amino_acid_list, i, sequence_length)
        if positive:
            bind_record = generate_binding_record(sequence_length, i, interval)
            pssm = generate_pssm(i, sequence_length, interval, random_flag)
            secondary_structure = generate_secondary_structure(bind_record, sequence_length)
        else:
            bind_record = "0" * sequence_length
            pssm = [[0]*20] * sequence_length
            secondary_structure = generate_secondary_structure(bind_record, sequence_length)

        write_to_bindres_file(bindres_dir, bindres_filename, proteinid, sequence, bind_record)
        write_to_pssm_file(pssm_dir, pssm_filename, pssm, sequence)
        write_to_secondary_structure_file(secondary_structure_dir, secondary_structure_filename, proteinid, secondary_structure, sequence)

        if positive:
            positive_or_negative = "positive"
        else:
            positive_or_negative = "negative"
        data_paths += "{}\t{}\t{}\t{}\n".format(proteinid, positive_or_negative, "bindres", "{}/{}".format(bindres_dir, bindres_filename))
        data_paths += "{}\t{}\t{}\t{}\n".format(proteinid, positive_or_negative, "pssm", "{}/{}".format(pssm_dir, pssm_filename))
        data_paths += "{}\t{}\t{}\t{}\n".format(proteinid, positive_or_negative, "secondary_structure", "{}/{}".format(secondary_structure_dir, secondary_structure_filename))
    return data_paths


def write_to_data_path_file(data_path_file, data_paths):
    with open(data_path_file, 'w') as fp:
        fp.write("#proteinid\tpositive_or_negative\tdatatype\tpath\n{}".format(data_paths))


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
        print Usage
        sys.exit(0)
    interval = int(sys.argv[1])
    if not interval > 0:
        raise ValueError("<interval> must be bigger than 0")
    sequence_length = int(sys.argv[2])
    if not sequence_length > 0:
        raise ValueError("<sequence_length> must be bigger than 0")
    random_flag = False
    if len(sys.argv) == 4 and sys.argv[3] == "--randomscore":
        random_flag = True
    sequence_length = int(sys.argv[2])

    workdir = "/Users/clclcocoro/work/lipid_bindResPred/program/test"
    bindres_dir = "{}/bindres".format(workdir)
    pssm_dir = "{}/pssm".format(workdir)
    secondary_structure_dir = "{}/secondary_structure".format(workdir)
    data_path_file = "{}/data_path_file.txt".format(workdir)
    # the same order in pssm
    amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    number_of_positive_proteins = 2
    number_of_negative_proteins = 3

    data_paths = generate_proteins_data(number_of_positive_proteins, bindres_dir, pssm_dir, secondary_structure_dir,
                                                                         sequence_length, interval, random_flag, positive=True)
    data_paths += generate_proteins_data(number_of_negative_proteins, bindres_dir, pssm_dir, secondary_structure_dir,
                                                                         sequence_length, interval, random_flag, positive=False)
    write_to_data_path_file(data_path_file, data_paths)
