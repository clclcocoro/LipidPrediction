#!/usr/bin/env python

"""Run binding prediction.


Usage:
  main.py -d <data_path_file> -w <window_size> -s <smoothing_window_size> -m <pickled_model_file> -o <output_file> [--original] [--exp] [--smoothed] [--exp_smoothed] [--aaindex] [--secondary_structure]
  main.py (-h | --help)
  main.py --version

Options:
  -d                           file that includes all paths of data
  -m                           pickled model file.
  -o                           output file.
  -w                           window size.
  -s                           smoothing window size.
  --original                   Use original pssm as feature.
  --exp                        Use exp pssm as feature.
  --smoothed                   Use smoothed pssm as feature.
  --exp_smoothed               Use exp smoothed pssm as feature.
  --aaindex                    Use AAindex as feature.
  --secondary_structure        Use secondary structure as feature.
  -h --help                    Show this screen.
  --version                    Show version.


data_path_file FORMAT (TSV)
#proteinid\tpositive_or_negative\tdatatype\tpath
proteinid1\tpositive\tbindres\t/Users/username/data/AAAA.bindres
proteinid1\tpositive\tpssm\t/Users/username/data/AAAA.pssm
proteinid1\tpositive\tsecondary_structure\t/Users/username/data/AAAA.secondary_structure.txt
proteinid2\tnegative\tbindres\t/Users/username/data/BBBB.bindres
proteinid2\tnegative\tpssm\t/Users/username/data/BBBB.pssm
proteinid2\tnegative\tsecondary_structure\t/Users/username/data/BBBB.secondary_structure.txt
...
...

The order of rows doesn't matter.



"""

from docopt import docopt
from sklearn import svm
import numpy
import filepath
import feature
import dataset
import pickle

def write_to_output_file(output_file, decision_values, predicted_labels, proteinid_index_list):
    with open(output_file, 'w') as fp:
        for proteinid_index, decision_value, predicted_label in zip(proteinid_index_list, decision_values, predicted_labels):
            fp.write('{}\t{}\t{}\t{}\n'.format(proteinid_index[0], proteinid_index[1], decision_value, predicted_label))


def create_ProteinHolder(smoothing_window_size, data_filepath):
    protein_holder = dataset.ProteinHolder()
    for proteinid in data_filepath.positive_proteinids:
        filepaths = data_filepath.get_filepaths_of_protein(proteinid)
        protein_holder.add_positive_protein(feature.Protein(proteinid, filepaths['pssm'], filepaths['secondary_structure'], filepaths['bindres'], smoothing_window_size=smoothing_window_size))
    for proteinid in data_filepath.negative_proteinids:
        filepaths = data_filepath.get_filepaths_of_protein(proteinid)
        protein_holder.add_negative_protein(feature.Protein(proteinid, filepaths['pssm'], filepaths['secondary_structure'], filepaths['bindres'], smoothing_window_size=smoothing_window_size))
    return protein_holder
 

if __name__ == "__main__":
    arguments = docopt(__doc__)
    data_path_file            = arguments['<data_path_file>']
    window_size               = int(arguments['<window_size>'])
    smoothing_window_size     = int(arguments['<smoothing_window_size>'])
    pickled_model_file        = arguments['<pickled_model_file>']
    output_file               = arguments['<output_file>']
    original_pssm             = arguments['--original']
    exp_pssm                  = arguments['--exp']
    smoothed_pssm             = arguments['--smoothed']
    exp_smoothed_pssm         = arguments['--exp_smoothed']
    AAindex                   = arguments['--aaindex']
    secondary_structure       = arguments['--secondary_structure']
    undersampling             = False
    shuffle                   = True

    # Normalizing AAindex.
    feature.normalize_AAindex()

    data_filepath = filepath.DataFilepath(data_path_file)
    protein_holder = create_ProteinHolder(smoothing_window_size, data_filepath)
    test_dataset, proteinid_index_list = dataset.create_test_dataset(protein_holder, window_size, original_pssm=original_pssm,
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm,
                                                       AAindex=AAindex, secondary_structure=secondary_structure, )
    pkl_file = open(pickled_model_file, 'rb')
    clf = pickle.load(pkl_file)
    pkl_file.close()
    decision_values = clf.decision_function(test_dataset)
    if type(decision_values[0]) is list or type(decision_values[0]) is numpy.ndarray:
        decision_values = map(lambda x: x[0], decision_values)
    predicted_labels = clf.predict(test_dataset)

    write_to_output_file(output_file, decision_values, predicted_labels, proteinid_index_list)
