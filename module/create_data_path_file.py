#!/usr/bin/env python

import os

if __name__ == "__main__":
    positive_pssm_dir                = "/Users/clclcocoro/work/lipid_bindResPred/work/pssm/positive"
    positive_bindres_dir             = "/Users/clclcocoro/work/lipid_bindResPred/work/bindres/positive"
    positive_secondary_structure_dir = "/Users/clclcocoro/work/lipid_bindResPred/jpred/secondary_struc/positive"
    negative_pssm_dir                = "/Users/clclcocoro/work/lipid_bindResPred/work/pssm/negative"
    negative_bindres_dir             = "/Users/clclcocoro/work/lipid_bindResPred/work/bindres/negative"
    negative_secondary_structure_dir = "/Users/clclcocoro/work/lipid_bindResPred/jpred/secondary_struc/negative"
    data_path_file = "/Users/clclcocoro/work/lipid_bindResPred/data_path_file.txt"
    
    with open(data_path_file, 'w') as fp:
        fp.write("#proteinid\tpositive_or_negative\tdatatype\tpath\n")
        for pssm_file in os.listdir(positive_pssm_dir):
            proteinid = pssm_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "positive", "pssm", "{}/{}".format(positive_pssm_dir, pssm_file)))
        for bindres_file in os.listdir(positive_bindres_dir):
            proteinid = bindres_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "positive", "bindres", "{}/{}".format(positive_bindres_dir, bindres_file)))
        for secondary_structure_file in os.listdir(positive_secondary_structure_dir):
            proteinid = secondary_structure_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "positive", "secondary_structure", "{}/{}".format(positive_secondary_structure_dir, secondary_structure_file)))
        for pssm_file in os.listdir(negative_pssm_dir):
            proteinid = pssm_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "negative", "pssm", "{}/{}".format(negative_pssm_dir, pssm_file)))
        for bindres_file in os.listdir(negative_bindres_dir):
            proteinid = bindres_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "negative", "bindres", "{}/{}".format(negative_bindres_dir, bindres_file)))
        for secondary_structure_file in os.listdir(negative_secondary_structure_dir):
            proteinid = secondary_structure_file[:6]
            fp.write("{}\t{}\t{}\t{}\n".format(proteinid, "negative", "secondary_structure", "{}/{}".format(negative_secondary_structure_dir, secondary_structure_file)))
