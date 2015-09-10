#!/usr/bin/env python

"""Run Lipid binding prediction with Cross Validation.


Usage:
  main.py -d <data_path_file> -w <window_size> -s <smoothing_window_size> -o <output_log_file> [--original] [--exp] [--smoothed] [--aaindex] [--secondary_structure]
  main.py (-h | --help)
  main.py --version

Options:
  -d                       file that includes all paths of data
  -o                       output log file.
  -w                       window size.
  -s                       smoothing window size.
  --original               Use original pssm as feature.
  --exp                    Use exp pssm as feature.
  --smoothed               Use smoothed pssm as feature.
  --aaindex                Use AAindex as feature.
  --secondary_structure    Use secondary structure as feature.
  -h --help                Show this screen.
  --version                Show version.


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
import filepath
import feature
import dataset
import cross_validation


def create_feature_info(window_size, smoothing_window_size, original, exp_pssm, smoothed_pssm, AAindex, secondary_structure):
    return "w:{}, sw:{}, original_pssm:{}, exp_pssm:{}, smoothed_pssm:{}, AAindex:{}, secondary_structure:{}".format(window_size, 
            smoothing_window_size, original, exp_pssm, smoothed_pssm, AAindex, secondary_structure)


def create_ProteinHolder(smoothing_window_size, data_filepath):
    protein_holder = dataset.ProteinHolder()
    for proteinid in data_filepath.positive_proteinids:
        filepaths = data_filepath.get_filepaths_of_protein(proteinid)
        protein_holder.add_positive_protein(feature.Protein(filepaths['pssm'], filepaths['secondary_structure'], filepaths['bindres']), smoothing_window_size=smoothing_window_size)
    for proteinid in data_filepath.negative_proteinids:
        filepaths = data_filepath.get_filepaths_of_protein(proteinid)
        protein_holder.add_negative_protein(feature.Protein(filepaths['pssm'], filepaths['secondary_structure'], filepaths['bindres']), smoothing_window_size=smoothing_window_size)
    return protein_holder
 

def grid_search(cost_range, gamma_range, crossValidation):
    for i in xrange(cost_range[0], cost_range[1]+1):
        for j in xrange(gamma_range[0], gamma_range[1]+1):
            cost  = 2 ** i
            gamma = 2 ** j
            cross_validation.run_SVM(cost, gamma)


if __name__ == "__main__"
    arguments = docopt(__doc__)
    data_path_file        = arguments['<data_path_file>']
    log_file              = arguments['<output_log_file>']
    window_size           = int(arguments['<window_size>'])
    smoothing_window_size = int(arguments['<smoothing_window_size>'])
    original_pssm         = arguments['--original']
    exp_pssm              = arguments['--exp']
    smoothed_pssm         = arguments['--smoothed']
    AAindex               = arguments['--aaindex']
    secondary_structure   = arguments['--secondary_structure']
    undersampling         = False
    shuffle               = True
    fold                  = 5
    cost_range            = (-4, 4) # 2**this_integer
    gamma_range           = (-4, 4) # 2**this_integer

    # Normalizing AAindex.
    feature.normalize_AAindex()

    data_filepath = filepath.DataFilepath(data_path_file)
    protein_holder = create_ProteinHolder(smoothing_window_size, data_filepath)
    feature_info = create_feature_info(window_size, smoothing_window_size, original, 
                                    exp_pssm, smoothed_pssm, AAindex, secondary_structure)
    positive_dataset = dataset.create_positive_dataset(protein_holder, window_size, original_pssm=original_pssm
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, 
                                                       AAindex=AAindex, secondary_structure=secondary_structure)
    negative_dataset = dataset.create_negative_dataset(protein_holder, window_size, original_pssm=original_pssm
                                                       exp_pssm=exp_pssm, smoothed_pssm=smoothed_pssm, 
                                                       AAindex=AAindex, secondary_structure=secondary_structure)
    folded_dataset = dataset.FoldedDataset(positive_dataset, negative_dataset, fold=fold, undersampling=undersampling, shuffle=shuffle)
    crossValidation = cross_validation.CrossValidation(folded_dataset, feature_info, log_file, fold=fold)

    grid_search(cost_range, gamma_range, crossValidation)
