#!/usr/bin/env python

import re

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

    def parse_pssm_file(self, pssm_file)
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


def create_pssm_feature_vectors(protein, window_size):
    """
        0 0 0 0 0 0 0 0 0 0 ... 0
        0 0 0 0 0 0 0 0 0 0 ... 0
        -1 -2 -1 3 ...         -3
    """
    feature_vectors = []
    seqlen = len(protein.pssm)
    for i in xrange(seqlen):
        feature_vector = []
        p = i - window_size 
        if p < 0:
            feature_vector += [0] * (20 * (window_size-i))
            p = 0
        if i + window_size <= seqlen - 1:
            while p <= i + window_size:
                feature_vector += protein.pssm[p]
                p += 1
        else:
            while p <= seqlen - 1:
                feature_vector += protein.pssm[p]
                p += 1
            feature_vector += [0] * (20*((i+window_size)-(seqlen-1)))
        feature_vectors.append(feature_vector)
    return feature_vectors


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
