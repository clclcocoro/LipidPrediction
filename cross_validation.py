#!/usr/bin/env python

from sklearn import svm
import numpy
import feature
import dataset
import validate_performance

class CrossValidation(object):

    """
    Gene Scale in GA has to be (0, n).  n is greater than 0.
    """
    def __init__(self, folded_dataset, feature_info, log_file, fold=5):
        self.folded_dataset = folded_dataset
        self.log_file = log_file
        self.fold = fold
        self.feature_info = feature_info

    def write_log(self, cost, gamma, mean_AUC, mean_decision_value, mean_mcc, mean_performance_by_decision_value, mean_performance_by_predict_function):
        with open(self.log_file, 'a') as fp:
            #fp.write("{} {} {} {}\n".format(gene1, gene2, gene3, mean_AUC))
            fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.feature_info, cost, gamma, mean_AUC, mean_decision_value, mean_mcc, '\t'.join(map(str, mean_performance_by_decision_value)), '\t'.join(map(str, mean_performance_by_predict_function))))

    def add_performances(self, mean, new):
        for i in xrange(len(mean)):
            mean[i] += new[i]
        return mean

    def run_SVM(self, cost, gamma):
        mean_AUC = 0
        mean_decision_value = 0
        mean_mcc = 0
        mean_performance_by_decision_value = [0, 0, 0, 0]   # [SE, SP, ACC, MCC]
        mean_performance_by_predict_function = [0, 0, 0, 0] # [SE, SP, ACC, MCC]
        for test_fold in xrange(self.fold):
            test_labels, test_dataset, train_labels, train_dataset = self.folded_dataset.get_test_and_training_dataset(test_fold)
            clf = svm.SVC(C=cost, gamma=gamma, class_weight='auto')
            clf.fit(train_dataset, train_labels)
            decision_values = clf.decision_function(test_dataset)
            if type(decision_values[0]) is list or type(decision_values[0]) is numpy.ndarray:
                decision_values = map(lambda x: x[0], decision_values)
            AUC, decision_value_and_max_mcc = validate_performance.calculate_AUC(decision_values, test_labels)
            predicted_labels_by_decision_value = [1 if decision_value >= decision_value_and_max_mcc[0] else 0 for decision_value in decision_values]
            predicted_labels_by_predict_function = clf.predict(test_dataset)
            # [SE, SP, ACC, MCC]
            performances = validate_performance.calculate_performance(test_labels, predicted_labels_by_decision_value)
            mean_performance_by_decision_value = self.add_performances(mean_performance_by_decision_value, performances)
            performances = validate_performance.calculate_performance(test_labels, predicted_labels_by_predict_function)
            mean_performance_by_predict_function = self.add_performances(mean_performance_by_predict_function, performances)
            mean_AUC += AUC
            mean_decision_value += decision_value_and_max_mcc[0]
            mean_mcc += decision_value_and_max_mcc[1]
        mean_AUC /= self.fold
        mean_decision_value /= self.fold
        mean_mcc /= self.fold
        mean_performance_by_decision_value = map(lambda x: x/float(self.fold), mean_performance_by_decision_value)
        mean_performance_by_predict_function = map(lambda x: x/float(self.fold), mean_performance_by_predict_function)
        self.write_log(cost, gamma, mean_AUC, mean_decision_value, mean_mcc,
                        mean_performance_by_decision_value, mean_performance_by_predict_function)
