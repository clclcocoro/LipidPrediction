#!/usr/bin/env python

import sys
sys.path.append("..")

import os
import math
import unittest
import validate_performance
import feature
import dataset
import filepath
import main

def create_positive_and_negative_dataset():
    os.system("python create_testdata.py 9 30")

# Create Test Data
create_positive_and_negative_dataset()
workdir = "/Users/clclcocoro/work/lipid_bindResPred/program/test"
data_path_file = "/Users/clclcocoro/work/lipid_bindResPred/program/test/data_path_file.txt"


class TestValidate(unittest.TestCase):

    def test_calculate_TPR_FPR(self):
        TP, FP, TN, FN = 1, 1, 1, 1
        FPR, TPR = validate_performance.calculate_TPR_FPR(TP, FP, TN, FN)
        self.assertEqual(TPR, 0.5)
        self.assertEqual(FPR, 0.5)

        TP, FP, TN, FN = 1, 2, 3, 4
        FPR, TPR = validate_performance.calculate_TPR_FPR(TP, FP, TN, FN)
        self.assertEqual(TPR, 0.2)
        self.assertEqual(FPR, 0.4)

    def test_calculate_AUC(self):
        decision_values = [-1, -0.5, -0.1,  0.1, 0.5, 1]
        correct_labels = [0, 0, 0, 1, 1, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 1.0)

        decision_values = [-1, 0.5, -0.1,  0.1, -0.5, 1]
        correct_labels = [0, 1, 0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 1.0)

        decision_values = [1.0]*6
        correct_labels = [0, 1, 0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 0.5)

        decision_values = [-1, -0.5, 0.5, 1]
        correct_labels = [0, 1, 0, 1]
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(AUC, 0.75)

        decision_values = [1.0] * 2 + [-1.0] * 8 + [1.0] * 2 + [0.5] * 8
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.82)

        decision_values = [-1.0] * 8 + [0.0] * 2 + [0.0] * 2 + [0.5] * 8
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.98)

        decision_values = [-1.0] * 8 + [0.2] * 2 + [0.1]*3 +[0.2]*2 + [0.5]*5
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.92)

        decision_values = [-1.0]*6 + [0.0]*2+ [0.2]*2 + [0.0]*2+ [0.1]*1 +[0.2]*2 + [0.5]*5
        correct_labels = [0] * 10 + [1] * 10
        AUC, mcc = validate_performance.calculate_AUC(decision_values, correct_labels)
        self.assertEqual(round(AUC*(10**5))/(10**5), 0.9)
 

class TestFilepath(unittest.TestCase):

    def test_data_filepath(self):
        c_positive_proteinids = ["PPPP:{}".format(i) for i in xrange(2)]
        c_positive_bindres_paths = map(lambda x: "{}/bindres/{}.bindres".format(workdir, x), c_positive_proteinids)
        c_positive_pssm_paths = map(lambda x: "{}/pssm/{}.pssm".format(workdir, x), c_positive_proteinids)
        c_positive_secondary_structure_paths = map(lambda x: "{}/secondary_structure/{}.secondary_structure.txt".format(workdir, x), c_positive_proteinids)
        c_negative_proteinids = ["NNNN:{}".format(i) for i in xrange(3)]
        c_negative_bindres_paths = map(lambda x: "{}/bindres/{}.bindres".format(workdir, x), c_negative_proteinids)
        c_negative_pssm_paths = map(lambda x: "{}/pssm/{}.pssm".format(workdir, x), c_negative_proteinids)
        c_negative_secondary_structure_paths = map(lambda x: "{}/secondary_structure/{}.secondary_structure.txt".format(workdir, x), c_negative_proteinids)

        data_filepath = filepath.DataFilepath(data_path_file)

        positive_proteinids = data_filepath.get_positive_proteinids()
        for c, t in zip(c_positive_proteinids, positive_proteinids):
            self.assertEqual(c, t)
        for i, positive_proteinid in enumerate(positive_proteinids):
            filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinid)
            self.assertEqual(c_positive_bindres_paths[i], filepath_of_protein['bindres'])
            self.assertEqual(c_positive_pssm_paths[i], filepath_of_protein['pssm'])
            self.assertEqual(c_positive_secondary_structure_paths[i], filepath_of_protein['secondary_structure'])

        negative_proteinids = data_filepath.get_negative_proteinids()
        for c, t in zip(c_negative_proteinids, negative_proteinids):
            self.assertEqual(c, t)
        for i, negative_proteinid in enumerate(negative_proteinids):
            filepath_of_protein = data_filepath.get_filepaths_of_protein(negative_proteinid)
            self.assertEqual(c_negative_bindres_paths[i], filepath_of_protein['bindres'])
            self.assertEqual(c_negative_pssm_paths[i], filepath_of_protein['pssm'])
            self.assertEqual(c_negative_secondary_structure_paths[i], filepath_of_protein['secondary_structure'])


class TestFeature(unittest.TestCase):

    def test_protein_all_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([1]*20+[-1]*20+[-1]*20)
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20+[-1]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.all_pssm_feature_vectors(window_size)
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        window_size = 2
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_pssm_feature_vectors.append([0]*20*window_size+[2]*20+[1]*20+[-1]*20)
            elif i == 1:
                c_pssm_feature_vectors.append([0]*20*1+[2]*20+[1]*20+[-1]*20*window_size)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20*window_size+[1]*20+[0]*20*window_size)
            elif i == sequence_length - 2:
                c_pssm_feature_vectors.append([-1]*20*(window_size+1)+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20+[-1]*20*window_size)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20*(window_size+1))
            elif i % (interval+1) == 3:
                c_pssm_feature_vectors.append([1]*20+[-1]*20*(window_size+2))
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20*window_size+[1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20*window_size+[-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_pssm_feature_vectors.append([-1]*20*(window_size+2)+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.all_pssm_feature_vectors(window_size)
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_bind_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) != 0:
                continue
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.bind_pssm_feature_vectors(window_size)
        self.assertEqual(3, len(pssm_feature_vectors))
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_non_bind_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) == 0:
                continue
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([1]*20+[-1]*20+[-1]*20)
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20+[-1]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
 
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.non_bind_pssm_feature_vectors(window_size)
        self.assertEqual(27, len(pssm_feature_vectors))
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_protein_all_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.all_smoothed_pssm_feature_vectors(window_size)
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        smoothing_window_size = 2
        window_size = 2
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20*window_size+[2]*20+[1]*20+[0]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20+[0]*20+[-3]*20)
            elif i == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[1]*20+[0]*20+[-3]*20+[-5]*20)
            elif i == 3:
                c_smoothed_pssm_feature_vectors.append([1]*20+[0]*20+[-3]*20+[-5]*20*2)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-2]*20+[-1]*20+[0]*20*window_size)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-5]*20+[-3]*20+[-2]*20+[-1]*20+[0]*20)
            elif i == sequence_length - 3:
                c_smoothed_pssm_feature_vectors.append([-5]*20*2+[-3]*20+[-2]*20+[-1]*20)
            elif i == sequence_length - 4:
                c_smoothed_pssm_feature_vectors.append([-5]*20*3+[-3]*20+[-2]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[2]*20+[2]*20+[2]*20+[0]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([2]*20+[2]*20+[2]*20+[0]*20+[-3]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[2]*20+[0]*20+[-3]*20+[-5]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([2]*20+[0]*20+[-3]*20+[-5]*20*window_size)
            elif i % (interval+1) == 4:
                c_smoothed_pssm_feature_vectors.append([0]*20+[-3]*20+[-5]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[0]*20+[2]*20+[2]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-5]*20+[-3]*20+[0]*20+[2]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-5]*20*window_size+[-3]*20+[0]*20+[2]*20)
            elif i % (interval+1) == interval-3:
                c_smoothed_pssm_feature_vectors.append([-5]*20*(window_size+1)+[-3]*20+[0])
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-5]*20*(window_size+1)+[-3]*20)
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.all_smoothed_pssm_feature_vectors(window_size)
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_bind_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) != 0:
                continue
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.bind_smoothed_pssm_feature_vectors(window_size)
        self.assertEqual(3, len(smoothed_pssm_feature_vectors))
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_non_bind_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) == 0:
                continue
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.non_bind_smoothed_pssm_feature_vectors(window_size)
        self.assertEqual(27, len(smoothed_pssm_feature_vectors))
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_protein_all_exp_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([1]*20+[-1]*20+[-1]*20)
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20+[-1]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
        c_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.all_exp_pssm_feature_vectors(window_size)
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        window_size = 2
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_pssm_feature_vectors.append([0]*20*window_size+[2]*20+[1]*20+[-1]*20)
            elif i == 1:
                c_pssm_feature_vectors.append([0]*20*1+[2]*20+[1]*20+[-1]*20*window_size)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20*window_size+[1]*20+[0]*20*window_size)
            elif i == sequence_length - 2:
                c_pssm_feature_vectors.append([-1]*20*(window_size+1)+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20+[-1]*20*window_size)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20*(window_size+1))
            elif i % (interval+1) == 3:
                c_pssm_feature_vectors.append([1]*20+[-1]*20*(window_size+2))
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20*window_size+[1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20*window_size+[-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_pssm_feature_vectors.append([-1]*20*(window_size+2)+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
        c_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.all_exp_pssm_feature_vectors(window_size)
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_bind_exp_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) != 0:
                continue
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
        c_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.bind_exp_pssm_feature_vectors(window_size)
        self.assertEqual(3, len(pssm_feature_vectors))
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_non_bind_exp_pssm_feature_vectors(self):
        smoothing_window_size = 3
        window_size = 1
        interval = 9
        sequence_length = 30
        c_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) == 0:
                continue
            if i == 0:
                c_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20)
            elif i == sequence_length - 1:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_pssm_feature_vectors.append([1]*20+[2]*20+[1]*20)
            elif i % (interval+1) == 1:
                c_pssm_feature_vectors.append([2]*20+[1]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_pssm_feature_vectors.append([1]*20+[-1]*20+[-1]*20)
            elif i % (interval+1) == interval:
                c_pssm_feature_vectors.append([-1]*20+[1]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_pssm_feature_vectors.append([-1]*20+[-1]*20+[1]*20)
            else:
                c_pssm_feature_vectors.append([-1]*20*(2*window_size+1))
        c_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_pssm_feature_vectors]
 
        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        pssm_feature_vectors = protein.non_bind_exp_pssm_feature_vectors(window_size)
        self.assertEqual(27, len(pssm_feature_vectors))
        for c, t in zip(c_pssm_feature_vectors, pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_protein_all_exp_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))
        c_smoothed_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_smoothed_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.all_exp_smoothed_pssm_feature_vectors(window_size)
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        smoothing_window_size = 2
        window_size = 2
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20*window_size+[2]*20+[1]*20+[0]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([0]*20+[2]*20+[1]*20+[0]*20+[-3]*20)
            elif i == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[1]*20+[0]*20+[-3]*20+[-5]*20)
            elif i == 3:
                c_smoothed_pssm_feature_vectors.append([1]*20+[0]*20+[-3]*20+[-5]*20*2)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-2]*20+[-1]*20+[0]*20*window_size)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-5]*20+[-3]*20+[-2]*20+[-1]*20+[0]*20)
            elif i == sequence_length - 3:
                c_smoothed_pssm_feature_vectors.append([-5]*20*2+[-3]*20+[-2]*20+[-1]*20)
            elif i == sequence_length - 4:
                c_smoothed_pssm_feature_vectors.append([-5]*20*3+[-3]*20+[-2]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[2]*20+[2]*20+[2]*20+[0]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([2]*20+[2]*20+[2]*20+[0]*20+[-3]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[2]*20+[0]*20+[-3]*20+[-5]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([2]*20+[0]*20+[-3]*20+[-5]*20*window_size)
            elif i % (interval+1) == 4:
                c_smoothed_pssm_feature_vectors.append([0]*20+[-3]*20+[-5]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[0]*20+[2]*20+[2]*20+[2]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-5]*20+[-3]*20+[0]*20+[2]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-5]*20*window_size+[-3]*20+[0]*20+[2]*20)
            elif i % (interval+1) == interval-3:
                c_smoothed_pssm_feature_vectors.append([-5]*20*(window_size+1)+[-3]*20+[0])
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-5]*20*(window_size+1)+[-3]*20)
        c_smoothed_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_smoothed_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.all_exp_smoothed_pssm_feature_vectors(window_size)
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_bind_exp_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) != 0:
                continue
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))
        c_smoothed_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_smoothed_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.bind_exp_smoothed_pssm_feature_vectors(window_size)
        self.assertEqual(3, len(smoothed_pssm_feature_vectors))
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_feature_non_bind_exp_smoothed_pssm_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence_length = 30
        c_smoothed_pssm_feature_vectors = []
        for i in xrange(sequence_length):
            if i % (interval+1) == 0:
                continue
            if i == 0:
                c_smoothed_pssm_feature_vectors.append([0]*20+[3]*20+[2]*20)
            elif i == 1:
                c_smoothed_pssm_feature_vectors.append([3]*20+[2]*20+[-1]*20)
            elif i == sequence_length - 1:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[0]*20+[0]*20)
            elif i == sequence_length - 2:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[0]*20)
            elif i % (interval+1) == 0:
                c_smoothed_pssm_feature_vectors.append([2]*20+[4]*20+[2]*20)
            elif i % (interval+1) == 1:
                c_smoothed_pssm_feature_vectors.append([4]*20+[2]*20+[-1]*20)
            elif i % (interval+1) == 2:
                c_smoothed_pssm_feature_vectors.append([2]*20+[-1]*20+[-3]*20)
            elif i % (interval+1) == 3:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[-3]*20*(window_size+1))
            elif i % (interval+1) == interval:
                c_smoothed_pssm_feature_vectors.append([-1]*20+[2]*20+[4]*20)
            elif i % (interval+1) == interval-1:
                c_smoothed_pssm_feature_vectors.append([-3]*20+[-1]*20+[2]*20)
            elif i % (interval+1) == interval-2:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(window_size+1)+[-1]*20)
            else:
                c_smoothed_pssm_feature_vectors.append([-3]*20*(2*window_size+1))
        c_smoothed_pssm_feature_vectors = [map(lambda x: 1/(1+math.exp(-x)), feature_vector) for feature_vector in c_smoothed_pssm_feature_vectors]

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        smoothed_pssm_feature_vectors = protein.non_bind_exp_smoothed_pssm_feature_vectors(window_size)
        self.assertEqual(27, len(smoothed_pssm_feature_vectors))
        for c, t in zip(c_smoothed_pssm_feature_vectors, smoothed_pssm_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_all_AAindex_feature_vectors(self):
        # Normalize AAindex
        feature.normalize_AAindex()
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence = "A"*10+"R"*10+"N"*10
        sequence_length = 30
        c_AAindex_feature_vectors = []
        def calculate_AAindex_feature(AAindexes, aa):
            feat = []
            for AAindex in AAindexes:
                feat.append(AAindex[aa])
            return feat
        for i in xrange(sequence_length):
            buff = []
            if i == 0:
                buff += [0.5]*len(feature.AAindexes)
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            elif i == sequence_length-1:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += [0.5]*len(feature.AAindexes)
            else:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            c_AAindex_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        AAindex_feature_vectors = protein.all_AAindex_feature_vectors(window_size)
        self.assertEqual(30, len(AAindex_feature_vectors))
        for c, t in zip(c_AAindex_feature_vectors, AAindex_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        smoothing_window_size = 1
        window_size = 2
        interval = 9
        sequence = "A"*10+"R"*10+"N"*10
        sequence_length = 30
        c_AAindex_feature_vectors = []
        def calculate_AAindex_feature(AAindexes, aa):
            feat = []
            for AAindex in AAindexes:
                feat.append(AAindex[aa])
            return feat
        for i in xrange(sequence_length):
            buff = []
            if i == 0:
                buff += [0.5]*len(feature.AAindexes)*2
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+2])
            elif i == 1:
                buff += [0.5]*len(feature.AAindexes)
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+2])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+3])
            elif i == sequence_length-1:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-2])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += [0.5]*len(feature.AAindexes)*2
            elif i == sequence_length-2:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-3])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-2])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += [0.5]*len(feature.AAindexes)
            else:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-2])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+2])
            c_AAindex_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        AAindex_feature_vectors = protein.all_AAindex_feature_vectors(window_size)
        self.assertEqual(30, len(AAindex_feature_vectors))
        for c, t in zip(c_AAindex_feature_vectors, AAindex_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)



    def test_bind_AAindex_feature_vectors(self):
        # Normalize AAindex
        feature.normalize_AAindex()
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence = "A"*10+"R"*10+"N"*10
        sequence_length = 30
        c_AAindex_feature_vectors = []
        def calculate_AAindex_feature(AAindexes, aa):
            feat = []
            for AAindex in AAindexes:
                feat.append(AAindex[aa])
            return feat
        for i in xrange(sequence_length):
            buff = []
            if i % (interval+1) != 0:
                continue
            if i == 0:
                buff += [0.5]*len(feature.AAindexes)
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            elif i == sequence_length-1:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += [0.5]*len(feature.AAindexes)
            else:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            c_AAindex_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        AAindex_feature_vectors = protein.bind_AAindex_feature_vectors(window_size)
        self.assertEqual(3, len(AAindex_feature_vectors))
        for c, t in zip(c_AAindex_feature_vectors, AAindex_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_non_bind_AAindex_feature_vectors(self):
        # Normalize AAindex
        feature.normalize_AAindex()
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        sequence = "A"*10+"R"*10+"N"*10
        sequence_length = 30
        c_AAindex_feature_vectors = []
        def calculate_AAindex_feature(AAindexes, aa):
            feat = []
            for AAindex in AAindexes:
                feat.append(AAindex[aa])
            return feat
        for i in xrange(sequence_length):
            buff = []
            if i % (interval+1) == 0:
                continue
            if i == 0:
                buff += [0.5]*len(feature.AAindexes)
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            elif i == sequence_length-1:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += [0.5]*len(feature.AAindexes)
            else:
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i-1])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i])
                buff += calculate_AAindex_feature(feature.AAindexes, sequence[i+1])
            c_AAindex_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        AAindex_feature_vectors = protein.non_bind_AAindex_feature_vectors(window_size)
        self.assertEqual(27, len(AAindex_feature_vectors))
        for c, t in zip(c_AAindex_feature_vectors, AAindex_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_all_secondary_structure_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        secondary_structure = "HE"+"-"*7+"EHE"+"-"*7+"EHE"+"-"*8
        sequence_length = 30
        c_secondary_structure_feature_vectors = []
        def calculate_secondary_structure_feature(ss):
            if ss == '-':
                return [1, 0, 0]
            elif ss == 'H':
                return [0, 1, 0]
            elif ss == 'E':
                return [0, 0, 1]
        for i in xrange(sequence_length):
            buff = []
            if i == 0:
                buff += [0]*3
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            elif i == sequence_length-1:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += [0]*3
            else:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            c_secondary_structure_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        secondary_structure_feature_vectors = protein.all_secondary_structure_feature_vectors(window_size)
        self.assertEqual(30, len(secondary_structure_feature_vectors))
        for c, t in zip(c_secondary_structure_feature_vectors, secondary_structure_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

        smoothing_window_size = 1
        window_size = 2
        interval = 9
        secondary_structure = "HE"+"-"*7+"EHE"+"-"*7+"EHE"+"-"*8
        sequence_length = 30
        c_secondary_structure_feature_vectors = []
        def calculate_secondary_structure_feature(ss):
            if ss == '-':
                return [1, 0, 0]
            elif ss == 'H':
                return [0, 1, 0]
            elif ss == 'E':
                return [0, 0, 1]
        for i in xrange(sequence_length):
            buff = []
            if i == 0:
                buff += [0]*3*2
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
                buff += calculate_secondary_structure_feature(secondary_structure[i+2])
            elif i == 1:
                buff += [0]*3
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
                buff += calculate_secondary_structure_feature(secondary_structure[i+2])
            elif i == sequence_length-1:
                buff += calculate_secondary_structure_feature(secondary_structure[i-2])
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += [0]*3*2
            elif i == sequence_length-2:
                buff += calculate_secondary_structure_feature(secondary_structure[i-2])
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
                buff += [0]*3
            else:
                buff += calculate_secondary_structure_feature(secondary_structure[i-2])
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
                buff += calculate_secondary_structure_feature(secondary_structure[i+2])
            c_secondary_structure_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        secondary_structure_feature_vectors = protein.all_secondary_structure_feature_vectors(window_size)
        self.assertEqual(30, len(secondary_structure_feature_vectors))
        for c, t in zip(c_secondary_structure_feature_vectors, secondary_structure_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_bind_secondary_structure_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        secondary_structure = "HE"+"-"*7+"EHE"+"-"*7+"EHE"+"-"*8
        sequence_length = 30
        c_secondary_structure_feature_vectors = []
        def calculate_secondary_structure_feature(ss):
            if ss == '-':
                return [1, 0, 0]
            elif ss == 'H':
                return [0, 1, 0]
            elif ss == 'E':
                return [0, 0, 1]
        for i in xrange(sequence_length):
            buff = []
            if i % (interval+1) != 0:
                continue
            if i == 0:
                buff += [0]*3
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            elif i == sequence_length-1:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += [0]*3
            else:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            c_secondary_structure_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        secondary_structure_feature_vectors = protein.bind_secondary_structure_feature_vectors(window_size)
        self.assertEqual(3, len(secondary_structure_feature_vectors))
        for c, t in zip(c_secondary_structure_feature_vectors, secondary_structure_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)

    def test_non_bind_secondary_structure_feature_vectors(self):
        smoothing_window_size = 1
        window_size = 1
        interval = 9
        secondary_structure = "HE"+"-"*7+"EHE"+"-"*7+"EHE"+"-"*8
        sequence_length = 30
        c_secondary_structure_feature_vectors = []
        def calculate_secondary_structure_feature(ss):
            if ss == '-':
                return [1, 0, 0]
            elif ss == 'H':
                return [0, 1, 0]
            elif ss == 'E':
                return [0, 0, 1]
        for i in xrange(sequence_length):
            buff = []
            if i % (interval+1) == 0:
                continue
            if i == 0:
                buff += [0]*3
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            elif i == sequence_length-1:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += [0]*3
            else:
                buff += calculate_secondary_structure_feature(secondary_structure[i-1])
                buff += calculate_secondary_structure_feature(secondary_structure[i])
                buff += calculate_secondary_structure_feature(secondary_structure[i+1])
            c_secondary_structure_feature_vectors.append(buff)

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_proteinids = data_filepath.get_positive_proteinids()
        filepath_of_protein = data_filepath.get_filepaths_of_protein(positive_proteinids.pop())
        protein = feature.Protein(filepath_of_protein['pssm'], filepath_of_protein['secondary_structure'], filepath_of_protein['bindres'], smoothing_window_size=smoothing_window_size)
        secondary_structure_feature_vectors = protein.non_bind_secondary_structure_feature_vectors(window_size)
        self.assertEqual(27, len(secondary_structure_feature_vectors))
        for c, t in zip(c_secondary_structure_feature_vectors, secondary_structure_feature_vectors):
            for c_ele, t_ele in zip(c, t):
                self.assertEqual(c_ele, t_ele)


class TestDataset(unittest.TestCase):

    def test_concatnate_feature_vector_list(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        b = [[10, 11, 12], [13, 14, 15],[16, 17, 18]]
        c_result = [[1, 2, 3, 10, 11, 12], [4, 5, 6, 13, 14, 15], [7, 8, 9, 16, 17, 18]]
        result = dataset.concatnate_feature_vector_list(a, b)
        for c_vec, vec in zip(c_result, result):
            for c_ele, ele in zip(c_vec, vec):
                self.assertEqual(c_ele, ele)

    def test_create_positive_dataset(self):
        smoothing_window_size = 1
        window_size           = 1
        original_pssm         = True
        exp_pssm              = True
        smoothed_pssm         = True
        exp_smoothed_pssm     = True
        AAindex               = True
        secondary_structure   = True

        data_filepath = filepath.DataFilepath(data_path_file)
        protein_holder = main.create_ProteinHolder(smoothing_window_size, data_filepath)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(294, len(positive_dataset[0]))

        window_size           = 1
        original_pssm         = True
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(60, len(positive_dataset[0]))
        self.assertEqual(0, positive_dataset[0][0])
        self.assertEqual(2, positive_dataset[0][20])
        self.assertEqual(1, positive_dataset[0][40])

        window_size           = 2
        original_pssm         = True
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(100, len(positive_dataset[0]))
        self.assertEqual(0, positive_dataset[0][0])
        self.assertEqual(0, positive_dataset[0][20])
        self.assertEqual(2, positive_dataset[0][40])
        self.assertEqual(1, positive_dataset[0][60])
        self.assertEqual(-1, positive_dataset[0][80])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = True
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(60, len(positive_dataset[0]))
        self.assertEqual(0.5, positive_dataset[0][0])
        self.assertEqual(1/(1+math.exp(-2)), positive_dataset[0][20])
        self.assertEqual(1/(1+math.exp(-1)), positive_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = True
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(60, len(positive_dataset[0]))
        self.assertEqual(0, positive_dataset[0][0])
        self.assertEqual(3, positive_dataset[0][20])
        self.assertEqual(2, positive_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = True
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(60, len(positive_dataset[0]))
        self.assertEqual(0.5, positive_dataset[0][0])
        self.assertEqual(1/(1+math.exp(-3)), positive_dataset[0][20])
        self.assertEqual(1/(1+math.exp(-2)), positive_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = True
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(45, len(positive_dataset[0]))
        self.assertEqual(0.5, positive_dataset[0][0])
        self.assertEqual(feature.AAindexes[0]['A'], positive_dataset[0][15])
        self.assertEqual(feature.AAindexes[0]['R'], positive_dataset[1][15])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = True

        data_filepath = filepath.DataFilepath(data_path_file)
        positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(6, len(positive_dataset))
        self.assertEqual(9, len(positive_dataset[0]))
        self.assertEqual(0, positive_dataset[0][0])
        self.assertEqual(1, positive_dataset[0][4])
        self.assertEqual(1, positive_dataset[0][8])

    def test_create_negative_dataset(self):
        smoothing_window_size = 1
        window_size           = 1
        original_pssm         = True
        exp_pssm              = True
        smoothed_pssm         = True
        exp_smoothed_pssm     = True
        AAindex               = True
        secondary_structure   = True

        data_filepath = filepath.DataFilepath(data_path_file)
        protein_holder = main.create_ProteinHolder(smoothing_window_size, data_filepath)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(294, len(negative_dataset[0]))

        window_size           = 1
        original_pssm         = True
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(60, len(negative_dataset[0]))
        self.assertEqual(0, negative_dataset[0][0])
        self.assertEqual(0, negative_dataset[0][20])
        self.assertEqual(0, negative_dataset[0][40])

        window_size           = 2
        original_pssm         = True
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(100, len(negative_dataset[0]))
        self.assertEqual(0, negative_dataset[0][0])
        self.assertEqual(0, negative_dataset[0][20])
        self.assertEqual(0, negative_dataset[0][40])
        self.assertEqual(0, negative_dataset[0][60])
        self.assertEqual(0, negative_dataset[0][80])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = True
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(60, len(negative_dataset[0]))
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][0])
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][20])
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = True
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(60, len(negative_dataset[0]))
        self.assertEqual(0, negative_dataset[0][0])
        self.assertEqual(0, negative_dataset[0][20])
        self.assertEqual(0, negative_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = True
        AAindex               = False
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(60, len(negative_dataset[0]))
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][0])
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][20])
        self.assertEqual(1/(1+math.exp(0)), negative_dataset[0][40])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = True
        secondary_structure   = False

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(45, len(negative_dataset[0]))
        self.assertEqual(0.5, negative_dataset[0][0])
        self.assertEqual(feature.AAindexes[0]['A'], negative_dataset[0][15])
        self.assertEqual(feature.AAindexes[0]['A'], negative_dataset[1][15])
        self.assertEqual(feature.AAindexes[0]['R'], negative_dataset[10][15])

        window_size           = 1
        original_pssm         = False
        exp_pssm              = False
        smoothed_pssm         = False
        exp_smoothed_pssm     = False
        AAindex               = False
        secondary_structure   = True

        data_filepath = filepath.DataFilepath(data_path_file)
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm, 
                    exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm, AAindex=AAindex, secondary_structure=secondary_structure)
        self.assertEqual(90, len(negative_dataset))
        self.assertEqual(9, len(negative_dataset[0]))
        self.assertEqual(0, negative_dataset[0][0])
        self.assertEqual(1, negative_dataset[0][3])
        self.assertEqual(1, negative_dataset[0][6])












        
















if __name__ == "__main__":
    unittest.main()
