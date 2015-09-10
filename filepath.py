#!/usr/bin/env python

import re
import os


class DataFilepath(object):

    """
    self.filepaths = {'proteinid1': {'pssm': '/path/to/pssm1',
                                          'bindres': '/path/to/bindres1',
                                          'secondary_structure': '/path/to/secondary_structure1',
                                         },
                           'proteinid2': {'pssm': '/path/to/pssm2',
                                          'bindres': '/path/to/bindres2',
                                          'secondary_structure': '/path/to/secondary_structure2',
                                         },
                           ...
    
                           }

    """
    def __init__(self, data_path_file):
        self.data_path_file = data_path_file
        self.filepaths, self.positive_proteinids, self.negative_proteinids = self.parse_data_path_file(data_path_file)

    def parse_data_path_file(self, data_path_file):
        filepaths = {}
        positive_proteinids = []
        negative_proteinids = []
        with open(data_path_file) as fp:
            for line in fp:
                if re.match(r"^#", line):
                    continue
                proteinid, positive_or_negative, datatype, path = line.rstrip().split('\t')
                if not positive_or_negative in {'positive', 'negative'}:
                    raise ValueError("invalid positive_or_negative {}".format(positive_or_negative))
                if not datatype in {'pssm', 'bindres', 'secondary_structure'}:
                    raise ValueError("invalid datatype {}".format(datatype))
                if not proteinid in filepaths:
                    filepaths[proteinid] = {}
                filepaths[proteinid][datatype] = path
                if positive_or_negative == 'positive':
                    positive_proteinids.append(proteinid)
                elif positive_or_negative == 'negative':
                    negative_proteinids.append(proteinid)
        return filepaths, positive_proteinids, negative_proteinids

    def get_positive_proteinids(self):
        return self.positive_proteinids

    def get_negative_proteinids(self):
        return self.negative_proteinids
        
    def get_filepaths_of_protein(self, proteinid):
        return self.filepaths[proteinid]
