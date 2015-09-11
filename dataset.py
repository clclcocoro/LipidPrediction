#!/usr/bin/env python

import random
import copy


class ProteinHolder(object):

    def __init__(self):
        self.positive_proteins = []
        self.negative_proteins = []
        self.test_proteins = []

    def add_positive_protein(self, positive_protein):
        self.positive_proteins.append(positive_protein)

    def add_negative_protein(self, negative_protein):
        self.negative_proteins.append(negative_protein)

    def add_test_protein(self, test_protein):
        self.test_proteins.append(test_protein)

    def clear_positive_proteins(self):
        self.positive_proteins = []

    def clear_negative_proteins(self):
        self.negative_proteins = []

    def clear_test_proteins(self):
        self.test_proteins = []


def concatnate_feature_vector_list(a_vectors, b_vectors):
    if len(a_vectors) == 0:
        return b_vectors
    if len(b_vectors) == 0:
        return a_vectors
    concatnated_vectors = []
    for a_vector, b_vector in zip(a_vectors, b_vectors):
        concatnated_vectors.append(a_vector + b_vector)
    return concatnated_vectors


def create_positive_dataset(protein_holder, window_size, original_pssm=False, exp_pssm=True, smoothed_pssm=True, exp_smoothed_pssm=True, AAindex=True, secondary_structure=True):
    positive_dataset = []
    for protein in protein_holder.positive_proteins:
        feature_vectors = []
        if original_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_pssm_feature_vectors(window_size))
        if exp_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_exp_pssm_feature_vectors(window_size))
        if smoothed_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_smoothed_pssm_feature_vectors(window_size))
        if exp_smoothed_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_exp_smoothed_pssm_feature_vectors(window_size))
        if AAindex:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_AAindex_feature_vectors(window_size))
        if secondary_structure:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.bind_secondary_structure_feature_vectors(window_size))
        positive_dataset += feature_vectors

    return positive_dataset


def create_negative_dataset(protein_holder, window_size, original_pssm=False, exp_pssm=True, smoothed_pssm=True, exp_smoothed_pssm=True, AAindex=True, secondary_structure=True):
    negative_dataset = []
    for protein in protein_holder.negative_proteins:
        feature_vectors = []
        if original_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_pssm_feature_vectors(window_size))
        if exp_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_exp_pssm_feature_vectors(window_size))
        if smoothed_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_smoothed_pssm_feature_vectors(window_size))
        if exp_smoothed_pssm:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_exp_smoothed_pssm_feature_vectors(window_size))
        if AAindex:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_AAindex_feature_vectors(window_size))
        if secondary_structure:
            feature_vectors = concatnate_feature_vector_list(feature_vectors, protein.all_secondary_structure_feature_vectors(window_size))
        negative_dataset += feature_vectors
    return negative_dataset


class FoldedDataset(object):

    def __init__(self, positive_dataset, negative_dataset, fold=5, undersampling=True, shuffle=True):
        self.fold = fold
        positive_size = len(positive_dataset)
        negative_size = len(negative_dataset)
        if shuffle:
            my_seed=12345
            random.seed(my_seed)
            random.shuffle(positive_dataset)
            random.seed(my_seed)
            random.shuffle(negative_dataset)
        if undersampling:
            if negative_size > 2 * positive_size:
                negative_dataset = negative_dataset[:2*positive_size]
            if positive_size > 2 * negative_size:
                positive_dataset = positive_dataset[:2*negative_size]
        self.positive_size = len(positive_dataset)
        self.negative_size = len(negative_dataset)
        self.folded_positive_dataset = self.folding(self.positive_size, positive_dataset, fold)
        self.folded_negative_dataset = self.folding(self.negative_size, negative_dataset, fold)

    def folding(self, size, dataset, fold):
        p = size / fold
        r = size % fold
        if r == 0:
            folded_dataset = [dataset[p*i:p*(i+1)] for i in xrange(fold)]
        else:
            folded_dataset = [dataset[p*i:p*(i+1)] for i in xrange(fold)]
            for i in xrange(r):
                folded_dataset[i] += [dataset[fold*p+i]]
        return folded_dataset

    def get_folded_positive_dataset(self):
        return self.folded_positive_dataset

    def get_folded_negative_dataset(self):
        return self.folded_negative_dataset

    def get_test_and_training_dataset(self, test_fold, undersampling_training_dataset=True):
        if test_fold < 0 or self.fold <= test_fold:
            raise ValueError("test_fold [{}] must be between 0 and self.fold-1 [{}]".format(test_fold, self.fold-1))
        test_positive_size = len(self.folded_positive_dataset[test_fold])
        test_negative_size = len(self.folded_negative_dataset[test_fold])
        test_labels = [1] * test_positive_size + [0] * test_negative_size
        test_dataset = self.folded_positive_dataset[test_fold] + self.folded_negative_dataset[test_fold]
        positive_train_labels = [1] * (self.positive_size - test_positive_size)
        negative_train_labels = [0] * (self.negative_size - test_negative_size)
        positive_train_dataset = []
        negative_train_dataset = []
        for i in xrange(self.fold):
            if i != test_fold:
                positive_train_dataset += self.folded_positive_dataset[i]
                negative_train_dataset += self.folded_negative_dataset[i]
        if undersampling_training_dataset:
            if self.positive_size - test_positive_size > 1500:
                positive_train_labels = [1] * 1500
                positive_train_dataset = positive_train_dataset[:1500]
            if self.negative_size - test_negative_size > 1500:
                negative_train_labels = [0] * 1500
                negative_train_dataset = negative_train_dataset[:1500]
        train_labels = positive_train_labels + negative_train_labels
        train_dataset = positive_train_dataset + negative_train_dataset
        return test_labels, test_dataset, train_labels, train_dataset
