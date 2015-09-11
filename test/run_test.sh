#!/bin/bash

python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original --exp 
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original --exp --smoothed 
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original --exp --smoothed --exp_smoothed 
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original --exp --smoothed --exp_smoothed --aaindex 
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --original --exp --smoothed --exp_smoothed --aaindex --secondary_structure
python ../main.py -d ./data_path_file.txt -w 1 -s 1 -o ./test_output.txt --exp_smoothed --aaindex --secondary_structure
python ../main.py -d ./data_path_file.txt -w 1 -s 2 -o ./test_output.txt --exp_smoothed --aaindex --secondary_structure
python ../main.py -d ./data_path_file.txt -w 1 -s 3 -o ./test_output.txt --exp_smoothed --aaindex --secondary_structure
python ../main.py -d ./data_path_file.txt -w 1 -s 4 -o ./test_output.txt --exp_smoothed --aaindex --secondary_structure
python ../main.py -d ./data_path_file.txt -w 2 -s 1 -o ./test_output.txt --exp_smoothed --aaindex --secondary_structure
