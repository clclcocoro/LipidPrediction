#!/bin/bash

python ../main.py -d ./data_path_file_train.txt -w 1 -s 1 -o ./output_train.txt --exp_smoothed
python ../main.py -d ./data_path_file_train.txt -w 1 -s 2 -o ./output_train.txt --exp_smoothed
python ../main.py -d ./data_path_file_train.txt -w 1 -s 3 -o ./output_train.txt --exp_smoothed

COST=32
GAMMA=0.001953125
WINDOW=2
S_WINDOW=1
python ../create_model.py -d ./data_path_file_train.txt -w $WINDOW -s $S_WINDOW -o ./model.pickle -c $COST -g $GAMMA --exp
python ../predict.py -d ./data_path_file_test.txt -w $WINDOW -s $S_WINDOW -m ./model.pickle -o ./prediction.txt --exp
