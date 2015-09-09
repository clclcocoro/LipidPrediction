#!/usr/bin/env python

import re

# AAindexs
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
AAindexs = [
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


class Protein(object):

    """
    pssm = [[-1, -3, 0, ... 1], [-4, -5, 2, ..., -6], ..., [2, -1, 0, ..., -4]]
    secondary_structure = '-----HHHHH----BBB ... CC----'
    binding_record = '00000100000001100000000...000110'
    sequence = 'AAARAVAAAASRARRLPPPLPL...PRPALKKD'
    pdbid = '3K5H:A'
    """

    def __init__(self, pssm_file, secondary_structure_file, binding_residue_file):
        self.pssm = self.parse_pssm_file(pssm_file)
        self.secondary_structure = self.parse_secondary_structure_file(selsecondary_structure_file)
        self.binding_record, self.sequence, self.pdbid = self.parse_binding_residue_file(binding_residue_file)

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
                    li.append(int(line[idx:idx+3].strip()))
                pssm.append(li)

    def parse_secondary_structure_file(self, secondary_structure_file):
        with open(secondary_structure_file) as fp:
            for line in fp:
                if line[0] in {'-', 'H', 'B', 'C'}:
                    return line.rstrip()

    def parse_binding_residue_file(self, binding_residue_file):
        pdbid = ''
        sequence = ''
        binding_record = ''
        with open(binding_residue_file) as fp:
            for line in fp:
                if line[0] == '>':
                    pdbid = line[1:].rstrip()
                elif re.match(r'[A-Z]', line[0]): 
                    sequence = line.rstrip()
                elif line[0] in {'0', '1'}:
                    binding_record = line.rstrip()
        return binding_record, sequence, pdbid

    def is_binding_residue(self, position):
        if self.binding_record[position] == '1':
            return True
        elif self.binding_record[position] == '0':
            return False
        else:
            raise ValueError("binding_record is invalid {}", self.binding_record[position])

    def all_pssm_feature_vectors(self, window_size):
        feature_vectors = []
        seqlen = len(self.pssm)
        for i in xrange(seqlen):
            feature_vector = []
            p = i - window_size 
            if p < 0:
                feature_vector += [0] * (20 * (window_size-i))
                p = 0
            if i + window_size <= seqlen - 1:
                while p <= i + window_size:
                    feature_vector += self.pssm[p]
                    p += 1
            else:
                while p <= seqlen - 1:
                    feature_vector += self.pssm[p]
                    p += 1
                feature_vector += [0] * (20*((i+window_size)-(seqlen-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors

    def bind_pssm_feature_vectors(self, window_size):
        feature_vectors = []
        seqlen = len(self.pssm)
        for i in xrange(seqlen):
            if not is_binding_residue(i):
                continue
            feature_vector = []
            p = i - window_size 
            if p < 0:
                feature_vector += [0] * (20 * (window_size-i))
                p = 0
            if i + window_size <= seqlen - 1:
                while p <= i + window_size:
                    feature_vector += self.pssm[p]
                    p += 1
            else:
                while p <= seqlen - 1:
                    feature_vector += self.pssm[p]
                    p += 1
                feature_vector += [0] * (20*((i+window_size)-(seqlen-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors

    def non_bind_pssm_feature_vectors(self, window_size):
        feature_vectors = []
        seqlen = len(self.pssm)
        for i in xrange(seqlen):
            if is_binding_residue(i):
                continue
            feature_vector = []
            p = i - window_size 
            if p < 0:
                feature_vector += [0] * (20 * (window_size-i))
                p = 0
            if i + window_size <= seqlen - 1:
                while p <= i + window_size:
                    feature_vector += self.pssm[p]
                    p += 1
            else:
                while p <= seqlen - 1:
                    feature_vector += self.pssm[p]
                    p += 1
                feature_vector += [0] * (20*((i+window_size)-(seqlen-1)))
            feature_vectors.append(feature_vector)
        return feature_vectors

    def all_AAindex_feature_vectors(self):
        feature_vectors = []
        for aa in self.sequence:
            feature_vector = []
            for AAindex in AAindexs:
                feature_vector.append(AAindex[aa])
            feature_vectors.append(feature_vector)
        return feature_vectors

    def bind_AAindex_feature_vectors(self):
        feature_vectors = []
        for i, aa in enumerate(self.sequence):
            if not is_binding_residue(i):
                continue
            feature_vector = []
            for AAindex in AAindexs:
                feature_vector.append(AAindex[aa])
            feature_vectors.append(feature_vector)
        return feature_vectors

    def non_bind_AAindex_feature_vectors(self):
        feature_vectors = []
        for i, aa in enumerate(self.sequence):
            if is_binding_residue(i):
                continue
            feature_vector = []
            for AAindex in AAindexs:
                feature_vector.append(AAindex[aa])
            feature_vectors.append(feature_vector)
        return feature_vectors

    def all_secondary_structure_feature_vectors(self):
        # '-': [1, 0, 0, 0], 'H': [0, 1, 0, 0], 'B': [0, 0, 1, 0], 'C': [0, 0, 0, 1]
        feature_vectors = []
        for struc in self.secondary_structure:
            if struc == '-':
                feature_vectors.append([1, 0, 0, 0])
            elif struc == 'H':
                feature_vectors.append([0, 1, 0, 0])
            elif struc == 'B':
                feature_vectors.append([0, 0, 1, 0])
            elif struc == 'C':
                feature_vectors.append([0, 0, 0, 1])
            else:
                raise ValueError("secondary_structure is invalid {}", struc)
        return feature_vecotrs

    def bind_secondary_structure_feature_vectors(self):
        # '-': [1, 0, 0, 0], 'H': [0, 1, 0, 0], 'B': [0, 0, 1, 0], 'C': [0, 0, 0, 1]
        feature_vectors = []
        for i, struc in enumerate(self.secondary_structure):
            if not is_binding_residue(i):
                continue
            if struc == '-':
                feature_vectors.append([1, 0, 0, 0])
            elif struc == 'H':
                feature_vectors.append([0, 1, 0, 0])
            elif struc == 'B':
                feature_vectors.append([0, 0, 1, 0])
            elif struc == 'C':
                feature_vectors.append([0, 0, 0, 1])
            else:
                raise ValueError("secondary_structure is invalid {}", struc)
        return feature_vecotrs

    def non_bind_secondary_structure_feature_vectors(self):
        # '-': [1, 0, 0, 0], 'H': [0, 1, 0, 0], 'B': [0, 0, 1, 0], 'C': [0, 0, 0, 1]
        feature_vectors = []
        for i, struc in enumerate(self.secondary_structure):
            if is_binding_residue(i):
                continue
            if struc == '-':
                feature_vectors.append([1, 0, 0, 0])
            elif struc == 'H':
                feature_vectors.append([0, 1, 0, 0])
            elif struc == 'B':
                feature_vectors.append([0, 0, 1, 0])
            elif struc == 'C':
                feature_vectors.append([0, 0, 0, 1])
            else:
                raise ValueError("secondary_structure is invalid {}", struc)
        return feature_vecotors


def create_training_data(protein, window_size):
    positive_data = []
    negative_data = []
    create_pssm_feature_vectors(protein, window_size)
    for i, feature_vector in enumerate(feature_vectors):
        if i in bindRecord:
            positive_data.append(feature_vector)
        else:
            if i in negative_data_index_set:
                negative_data.append(feature_vector)
    return positive_data, negative_data


def create_dataset(bindingResidueData, pssmData, window_size):
    positive_dataset = []
    negative_dataset = []
    for uniprotURI in bindingResidueData.get_uniprotURIs():
        pssm = pssmData.get_PSSMRecord(uniprotURI)
        feature_vectors = create_feature_vectors(pssm, window_size)
        bindRecord = bindingResidueData.get_bindRecord(uniprotURI)
        negative_data_index_set = get_negative_data_index_set(bindRecord, len(pssm.get_PSSM()))
        positive_data, negative_data = create_training_data(bindRecord, feature_vectors, negative_data_index_set)
        positive_dataset += positive_data
        negative_dataset += negative_data
    return positive_dataset, negative_dataset


if __name__ == "__main__":
    bindres_file = "/Users/clclcocoro/galaxy/work/data/bindingData.txt"
    pssms_file = "/Users/clclcocoro/galaxy/work/data/pssms.txt"
    bindingResidueData, pssmData = parse_record_files(bindres_file, pssms_file)
    
    """
    print "bindingResidueData"
    print bindingResidueData
    print bindingResidueData.get_uniprotURIs()
    print bindingResidueData.get_bindRecord(bindingResidueData.get_uniprotURIs()[0])
    print "pssmData"
    print pssmData
    print pssmData.get_uniprotURIs()
    print pssmData.get_PSSMRecord(pssmData.get_uniprotURIs()[0]).get_PSSM()[:10]
    """
    
    window_size = 3
    positive_dataset, negative_dataset = create_dataset(bindingResidueData, pssmData, window_size)
    
    print "positive_dataset"
    print positive_dataset
    print "negative_dataset"
    print negative_dataset
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
