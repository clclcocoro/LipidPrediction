#!/usr/bin/env python

"""Run binding prediction with Cross Validation.


Usage:
  main.py -d <data_path_file> -w <window_size> -s <smoothing_window_size> -o <output_pickled_file> -c <cost> -g <gamma> [--original] [--exp] [--smoothed] [--exp_smoothed] [--aaindex] [--secondary_structure] [--only_positive_proteins] [--undersampling]
  main.py (-h | --help)
  main.py --version

Options:
  -d                           file that includes all paths of data
  -o                           output log file.
  -w                           window size.
  -s                           smoothing window size.
  --original                   Use original pssm as feature.
  --exp                        Use exp pssm as feature.
  --smoothed                   Use smoothed pssm as feature.
  --exp_smoothed               Use exp smoothed pssm as feature.
  --aaindex                    Use AAindex as feature.
  --secondary_structure        Use secondary structure as feature.
  --only_positive_proteins     Negative dataset extracted from binding proteins.
  --undersampling              Undersampling.
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

def create_SVM_classifier(cost, gamma, output_pickled_file, train_labels, train_dataset):
    clf = svm.SVC(C=cost, gamma=gamma, class_weight='auto')
    clf.fit(train_dataset, train_labels)
    with open(output_pickled_file, 'wb') as fp:
        pickle.dump(clf, fp)


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
    output_pickled_file       = arguments['<output_pickled_file>']
    cost                      = float(arguments['<cost>'])
    gamma                     = float(arguments['<gamma>'])
    original_pssm             = arguments['--original']
    exp_pssm                  = arguments['--exp']
    smoothed_pssm             = arguments['--smoothed']
    exp_smoothed_pssm         = arguments['--exp_smoothed']
    AAindex                   = arguments['--aaindex']
    secondary_structure       = arguments['--secondary_structure']
    only_positive_proteins    = arguments['--only_positive_proteins']
    undersampling             = arguments['--undersampling']
    shuffle                   = True
    fold                      = 5

    # Normalizing AAindex.
    feature.normalize_AAindex()

    data_filepath = filepath.DataFilepath(data_path_file)
    protein_holder = create_ProteinHolder(smoothing_window_size, data_filepath)
    positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm,
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm,
                                                       AAindex=AAindex, secondary_structure=secondary_structure)
    if only_positive_proteins:
        negative_dataset = dataset.create_negative_dataset_from_binding_proteins(protein_holder, window_size, original_pssm=original_pssm,
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm,
                                                       AAindex=AAindex, secondary_structure=secondary_structure)
    else:
        negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm,
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, exp_smoothed_pssm=exp_smoothed_pssm,
                                                       AAindex=AAindex, secondary_structure=secondary_structure)
    folded_dataset = dataset.FoldedDataset(positive_dataset, negative_dataset, fold=fold, undersampling=undersampling, shuffle=shuffle)
    train_labels, train_dataset = folded_dataset.get_training_dataset()
    create_SVM_classifier(cost, gamma, output_pickled_file, train_labels, train_dataset)
