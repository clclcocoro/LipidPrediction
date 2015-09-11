#!/usr/bin/env python

import re
import math

# AAindexes
# 
# H CIDH920101
# D Normalized hydrophobicity scales for alpha-proteins (Cid et al., 1992)
# H CIDH920102
# D Normalized hydrophobicity scales for beta-proteins (Cid et al., 1992)
# H CIDH920103
# D Normalized hydrophobicity scales for alpha+beta-proteins (Cid et al., 1992)
# H CIDH920104
# D Normalized hydrophobicity scales for alpha/beta-proteins (Cid et al., 1992)
# H CIDH920105
# D Normalized average hydrophobicity scales (Cid et al., 1992)
# H BHAR880101
# D Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
# H CHAM820101
# D Polarizability parameter (Charton-Charton, 1982)
# H CHAM820102
# D Free energy of solution in water, kcal/mole (Charton-Charton, 1982)
# H CHOC760101
# D Residue accessible surface area in tripeptide (Chothia, 1976)
# H CHOC760102
# D Residue accessible surface area in folded protein (Chothia, 1976)
# H BIGC670101
# D Residue volume (Bigelow, 1967)
# H CHAM810101
# D Steric parameter (Charton, 1981)
# H DAYM780101
# D Amino acid composition (Dayhoff et al., 1978a)
# H DAYM780201
# D Relative mutability (Dayhoff et al., 1978b)
# H FAUJ880102
# D Smoothed upsilon steric parameter (Fauchere et al., 1988)
AAindexes = [
    { 'A':  -0.450,'R':  -0.240,'N':  -0.200,'D':  -1.520,'C':   0.790,'Q':  -0.990,'E':  -0.800,'G':  -1.000,'H':   1.070,'I':   0.760,
      'L':   1.290,'K':  -0.360,'M':   1.370,'F':   1.480,'P':  -0.120,'S':  -0.980,'T':  -0.700,'W':   1.380,'Y':   1.490,'V':   1.260},
    { 'A':  -0.080,'R':  -0.090,'N':  -0.700,'D':  -0.710,'C':   0.760,'Q':  -0.400,'E':  -1.310,'G':  -0.840,'H':   0.430,'I':   1.390,
      'L':   1.240,'K':  -0.090,'M':   1.270,'F':   1.530,'P':  -0.010,'S':  -0.930,'T':  -0.590,'W':   2.250,'Y':   1.530,'V':   1.090},
    { 'A':   0.360,'R':  -0.520,'N':  -0.900,'D':  -1.090,'C':   0.700,'Q':  -1.050,'E':  -0.830,'G':  -0.820,'H':   0.160,'I':   2.170,
      'L':   1.180,'K':  -0.560,'M':   1.210,'F':   1.010,'P':  -0.060,'S':  -0.600,'T':  -1.200,'W':   1.310,'Y':   1.050,'V':   1.210},
    { 'A':   0.170,'R':  -0.700,'N':  -0.900,'D':  -1.050,'C':   1.240,'Q':  -1.200,'E':  -1.190,'G':  -0.570,'H':  -0.250,'I':   2.060,
      'L':   0.960,'K':  -0.620,'M':   0.600,'F':   1.290,'P':  -0.210,'S':  -0.830,'T':  -0.620,'W':   1.510,'Y':   0.660,'V':   1.210},
    { 'A':   0.020,'R':  -0.420,'N':  -0.770,'D':  -1.040,'C':   0.770,'Q':  -1.100,'E':  -1.140,'G':  -0.800,'H':   0.260,'I':   1.810,
      'L':   1.140,'K':  -0.410,'M':   1.000,'F':   1.350,'P':  -0.090,'S':  -0.970,'T':  -0.770,'W':   1.710,'Y':   1.110,'V':   1.130},
    { 'A':   0.357,'R':   0.529,'N':   0.463,'D':   0.511,'C':   0.346,'Q':   0.493,'E':   0.497,'G':   0.544,'H':   0.323,'I':   0.462,
      'L':   0.365,'K':   0.466,'M':   0.295,'F':   0.314,'P':   0.509,'S':   0.507,'T':   0.444,'W':   0.305,'Y':   0.420,'V':   0.386},
    { 'A':   0.046,'R':   0.291,'N':   0.134,'D':   0.105,'C':   0.128,'Q':   0.180,'E':   0.151,'G':   0.000,'H':   0.230,'I':   0.186,
      'L':   0.186,'K':   0.219,'M':   0.221,'F':   0.290,'P':   0.131,'S':   0.062,'T':   0.108,'W':   0.409,'Y':   0.298,'V':   0.140},
    { 'A':  -0.368,'R':  -1.030,'N':   0.000,'D':   2.060,'C':   4.530,'Q':   0.731,'E':   1.770,'G':  -0.525,'H':   0.000,'I':   0.791,
      'L':   1.070,'K':   0.000,'M':   0.656,'F':   1.060,'P':  -2.240,'S':  -0.524,'T':   0.000,'W':   1.600,'Y':   4.910,'V':   0.401},
    { 'A': 115.000,'R': 225.000,'N': 160.000,'D': 150.000,'C': 135.000,'Q': 180.000,'E': 190.000,'G':  75.000,'H': 195.000,'I': 175.000,
      'L': 170.000,'K': 200.000,'M': 185.000,'F': 210.000,'P': 145.000,'S': 115.000,'T': 140.000,'W': 255.000,'Y': 230.000,'V': 155.000},
    { 'A':  25.000,'R':  90.000,'N':  63.000,'D':  50.000,'C':  19.000,'Q':  71.000,'E':  49.000,'G':  23.000,'H':  43.000,'I':  18.000,
      'L':  23.000,'K':  97.000,'M':  31.000,'F':  24.000,'P':  50.000,'S':  44.000,'T':  47.000,'W':  32.000,'Y':  60.000,'V':  18.000},
    { 'A':  52.600,'R': 109.100,'N':  75.700,'D':  68.400,'C':  68.300,'Q':  89.700,'E':  84.700,'G':  36.300,'H':  91.900,'I': 102.000,
      'L': 102.000,'K': 105.100,'M':  97.700,'F': 113.900,'P':  73.600,'S':  54.900,'T':  71.200,'W': 135.400,'Y': 116.200,'V':  85.100},
    { 'A':   0.520,'R':   0.680,'N':   0.760,'D':   0.760,'C':   0.620,'Q':   0.680,'E':   0.680,'G':   0.000,'H':   0.700,'I':   1.020,
      'L':   0.980,'K':   0.680,'M':   0.780,'F':   0.700,'P':   0.360,'S':   0.530,'T':   0.500,'W':   0.700,'Y':   0.700,'V':   0.760},
    { 'A':   8.600,'R':   4.900,'N':   4.300,'D':   5.500,'C':   2.900,'Q':   3.900,'E':   6.000,'G':   8.400,'H':   2.000,'I':   4.500,
      'L':   7.400,'K':   6.600,'M':   1.700,'F':   3.600,'P':   5.200,'S':   7.000,'T':   6.100,'W':   1.300,'Y':   3.400,'V':   6.600},
    { 'A': 100.000,'R':  65.000,'N': 134.000,'D': 106.000,'C':  20.000,'Q':  93.000,'E': 102.000,'G':  49.000,'H':  66.000,'I':  96.000,
      'L':  40.000,'K':  56.000,'M':  94.000,'F':  41.000,'P':  56.000,'S': 120.000,'T':  97.000,'W':  18.000,'Y':  41.000,'V':  74.000},
    { 'A':   0.530,'R':   0.690,'N':   0.580,'D':   0.590,'C':   0.660,'Q':   0.710,'E':   0.720,'G':   0.000,'H':   0.640,'I':   0.960,
      'L':   0.920,'K':   0.780,'M':   0.770,'F':   0.710,'P':   0.000,'S':   0.550,'T':   0.630,'W':   0.840,'Y':   0.710,'V':   0.890},
    ]

# Replace above AAindexes to Normalize AAindexes
def normalize_AAindex():
    for i, AAindex in enumerate(AAindexes):
        keys = AAindex.keys()
        values = AAindex.values()
        max_val = max(values)
        min_val = min(values)
        normalized_values = map(lambda x: (x-min_val)/(max_val-min_val), values)
        for k, v in zip(keys, normalized_values):
            AAindex[k] = v
        AAindexes[i] = AAindex


def AAindex_feature(aa):
    feature = []
    for AAindex in AAindexes:
        feature.append(AAindex[aa])
    return feature
 

def secondary_structure_encode(struc):
    if struc == '-':
        return [1, 0, 0]
    elif struc == 'H':
        return [0, 1, 0]
    elif struc == 'E':
        return [0, 0, 1]
    else:
        raise ValueError("secondary_structure is invalid {}", struc)


class Protein(object):

    """
    pssm = [[-1, -3, 0, ... 1], [-4, -5, 2, ..., -6], ..., [2, -1, 0, ..., -4]]
    secondary_structure = '-----HHHHH----EEE ... HH----'
    binding_record = '00000100000001100000000...000110'
    sequence = 'AAARAVAAAASRARRLPPPLPL...PRPALKKD'
    proteinid = '3K5H:A'
    """

    def __init__(self, pssm_file, secondary_structure_file, binding_residue_file, smoothing_window_size=3):
        self.pssm = self.parse_pssm_file(pssm_file)
        self.secondary_structure = self.parse_secondary_structure_file(secondary_structure_file)
        self.binding_record, self.sequence, self.proteinid = self.parse_binding_residue_file(binding_residue_file)
        self.sequence_length = len(self.sequence)
        self.smoothing_window_size = smoothing_window_size
        self.smoothed_pssm = self.smoothe(smoothing_window_size)
        self.exp_pssm = [map(lambda x: 1/(1+math.exp(-x)), row) for row in self.pssm]
        self.exp_smoothed_pssm = [map(lambda x: 1/(1+math.exp(-x)), row) for row in self.smoothed_pssm]

    def parse_pssm_file(self, pssm_file):
        pssm = []
        with open(pssm_file) as fp:
            for c, line in enumerate(fp):
                if c <= 2: # Header
                    continue
                if not re.match(r'\s+\d+', line): # End of the Matrix.
                    return pssm
                li = []
                for i in xrange(9, 67, 3):
                    li.append(int(line[i:i+3].strip()))
                pssm.append(li)

    def parse_secondary_structure_file(self, secondary_structure_file):
        with open(secondary_structure_file) as fp:
            for c, line in enumerate(fp):
                if c == 2 and line[0] in {'-', 'H', 'E'}:
                    return line.rstrip()

    def parse_binding_residue_file(self, binding_residue_file):
        proteinid = ''
        sequence = ''
        binding_record = ''
        with open(binding_residue_file) as fp:
            for line in fp:
                if line[0] == '>':
                    proteinid = line[1:].rstrip()
                elif re.match(r'[A-Z]', line[0]): 
                    sequence = line.rstrip()
                elif line[0] in {'0', '1'}:
                    binding_record = line.rstrip()
        return binding_record, sequence, proteinid

    def init_smoothed_pssm(self, smoothing_window_size):
        self.smoothing_window_size = smoothing_window_size
        self.smoothed_pssm = self.smoothe(smoothing_window_size)

    def smoothe(self, smoothing_window_size):
        smoothed_pssm = []
        for i in xrange(self.sequence_length):
            smoothed_row = [0] * 20
            p = i - smoothing_window_size 
            if p < 0:
                p = 0
            if i + smoothing_window_size <= self.sequence_length - 1:
                while p <= i + smoothing_window_size:
                    for j in xrange(20):
                        smoothed_row[j] += self.pssm[p][j]
                    p += 1
            else:
                while p <= self.sequence_length - 1:
                    for j in xrange(20):
                        smoothed_row[j] += self.pssm[p][j]
                    p += 1
            smoothed_pssm.append(smoothed_row)
        return smoothed_pssm

    def is_binding_residue(self, position):
        if self.binding_record[position] == '1':
            return True
        elif self.binding_record[position] == '0':
            return False
        else:
            raise ValueError("binding_record is invalid {}", self.binding_record[position])

    def create_feature_vectors_from_pssm(self, pssm, window_size, exp_pssm=False, dataset_type='all'):
        if not dataset_type in {'all', 'bind', 'non_bind'}:
            raise ValueError("dataset_type is invalid {}".format(dataset_type))
        feature_vectors = []
        for i in xrange(self.sequence_length):
            if dataset_type == 'bind' and not self.is_binding_residue(i): 
                continue # Skip non-binding residue
            elif dataset_type == 'non_bind' and self.is_binding_residue(i):
                continue # Skip binding residue
            feature_vector = []
            p = i - window_size 
            if p < 0:
                if exp_pssm:
                    feature_vector += [0.5] * (20 * (window_size-i))
                else:
                    feature_vector += [0] * (20 * (window_size-i))
                p = 0
            if i + window_size <= self.sequence_length - 1:
                while p <= i + window_size:
                    feature_vector += pssm[p]
                    p += 1
            else:
                while p <= self.sequence_length - 1:
                    feature_vector += pssm[p]
                    p += 1
                if exp_pssm:
                    feature_vector += [0.5] * (20*((i+window_size)-(self.sequence_length-1)))
                else:
                    feature_vector += [0] * (20*((i+window_size)-(self.sequence_length-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors
        
    def all_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.pssm, window_size, dataset_type='all')

    def bind_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.pssm, window_size, dataset_type='bind')

    def non_bind_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.pssm, window_size, dataset_type='non_bind')
 
    def all_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.smoothed_pssm, window_size, dataset_type='all')

    def bind_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.smoothed_pssm, window_size, dataset_type='bind')

    def non_bind_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.smoothed_pssm, window_size, dataset_type='non_bind')

    def all_exp_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_pssm, window_size, exp_pssm=True, dataset_type='all')

    def bind_exp_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_pssm, window_size, exp_pssm=True, dataset_type='bind')

    def non_bind_exp_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_pssm, window_size, exp_pssm=True, dataset_type='non_bind')
 
    def all_exp_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_smoothed_pssm, window_size, exp_pssm=True, dataset_type='all')

    def bind_exp_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_smoothed_pssm, window_size, exp_pssm=True, dataset_type='bind')

    def non_bind_exp_smoothed_pssm_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_pssm(self.exp_smoothed_pssm, window_size, exp_pssm=True, dataset_type='non_bind')

    def create_feature_vectors_from_AAindex(self, window_size, dataset_type='all', normalized=True):
        if not dataset_type in {'all', 'bind', 'non_bind'}:
            raise ValueError("dataset_type is invalid {}".format(dataset_type))
        feature_vectors = []
        for i in xrange(self.sequence_length):
            if dataset_type == 'bind' and not self.is_binding_residue(i): 
                continue # Skip non-binding residue
            elif dataset_type == 'non_bind' and self.is_binding_residue(i):
                continue # Skip binding residue
            feature_vector = []
            p = i - window_size 
            if p < 0:
                if normalized:
                    feature_vector += [0.5] * (len(AAindexes) * (window_size-i))
                else:
                    feature_vector += [0] * (len(AAindexes) * (window_size-i))
                p = 0
            if i + window_size <= self.sequence_length - 1:
                while p <= i + window_size:
                    feature_vector += AAindex_feature(self.sequence[p])
                    p += 1
            else:
                while p <= self.sequence_length - 1:
                    feature_vector += AAindex_feature(self.sequence[p])
                    p += 1
                if normalized:
                    feature_vector += [0.5] * (len(AAindexes)*((i+window_size)-(self.sequence_length-1)))
                else:
                    feature_vector += [0] * (len(AAindexes)*((i+window_size)-(self.sequence_length-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors

    def all_AAindex_feature_vectors(self, window_size, normalized=True):
        return self.create_feature_vectors_from_AAindex(window_size, dataset_type='all', normalized=normalized)

    def bind_AAindex_feature_vectors(self, window_size, normalized=True):
        return self.create_feature_vectors_from_AAindex(window_size, dataset_type='bind', normalized=normalized)

    def non_bind_AAindex_feature_vectors(self, window_size, normalized=True):
        return self.create_feature_vectors_from_AAindex(window_size, dataset_type='non_bind', normalized=normalized)

    def create_feature_vectors_from_secondary_structure(self, window_size, dataset_type='all'):
        # '-': [1, 0, 0], 'H': [0, 1, 0], 'E': [0, 0, 1]
        if not dataset_type in {'all', 'bind', 'non_bind'}:
            raise ValueError("dataset_type is invalid {}".format(dataset_type))
        feature_vectors = []
        for i in xrange(self.sequence_length):
            if dataset_type == 'bind' and not self.is_binding_residue(i): 
                continue # Skip non-binding residue
            elif dataset_type == 'non_bind' and self.is_binding_residue(i):
                continue # Skip binding residue
            feature_vector = []
            p = i - window_size 
            if p < 0:
                feature_vector += [0] * (3 * (window_size-i))
                p = 0
            if i + window_size <= self.sequence_length - 1:
                while p <= i + window_size:
                    feature_vector += secondary_structure_encode(self.secondary_structure[p])
                    p += 1
            else:
                while p <= self.sequence_length - 1:
                    feature_vector += secondary_structure_encode(self.secondary_structure[p])
                    p += 1
                feature_vector += [0] * (3*((i+window_size)-(self.sequence_length-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors

    def all_secondary_structure_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_secondary_structure(window_size, dataset_type='all')

    def bind_secondary_structure_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_secondary_structure(window_size, dataset_type='bind')

    def non_bind_secondary_structure_feature_vectors(self, window_size):
        return self.create_feature_vectors_from_secondary_structure(window_size, dataset_type='non_bind')


if __name__ == "__main__":
    pssm_file = "/Users/clclcocoro/work/lipid_bindResPred/work/pssm/positive/1A25:B.pssm"
    bindres_file = "/Users/clclcocoro/work/lipid_bindResPred/work/bindres/positive/1A25:B.bindres"
    struc_file = "/Users/clclcocoro/work/lipid_bindResPred/jpred/secondary_struc/1A25:B.secondary_structure.txt"
    protein = Protein(pssm_file, struc_file, bindres_file)
    for row in protein.pssm:
        buff = ''
        for ele in row:
            buff += "{:3d}".format(ele)
        print buff
    print ""
    for row in protein.smoothed_pssm:
        buff = ''
        for ele in row:
            buff += "{:3d}".format(ele)
        print buff
    print len(protein.pssm)
    print len(protein.smoothed_pssm)
    print protein.secondary_structure
    print protein.binding_record
    print protein.proteinid
    print protein.sequence
    print len(protein.bind_pssm_feature_vectors(1))
    print protein.bind_AAindex_feature_vectors(1)
    print len(protein.bind_secondary_structure_feature_vectors(1))
    print len(protein.non_bind_pssm_feature_vectors(1))
    print len(protein.non_bind_AAindex_feature_vectors(1))
    print len(protein.non_bind_secondary_structure_feature_vectors(1))
    normalize_AAindex()
    print protein.bind_AAindex_feature_vectors(1)
    print len(protein.non_bind_AAindex_feature_vectors(1))
    #print len(protein.all_pssm_feature_vectors(1))
    #print len(protein.all_AAindex_feature_vectors(1))
    #print len(protein.all_secondary_structure_feature_vectors(1))
