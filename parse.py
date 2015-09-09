import re
import sys
import copy

aa = [
    'A', 
    'R',
    'N',
    'D',
    'C',
    'Q',
    'E',
    'G',
    'H',
    'I',
    'L',
    'K',
    'M',
    'F',
    'P',
    'S',
    'T',
    'W',
    'Y',
    'V'
    ]
flag = False
with open(sys.argv[1]) as fp:
    score = []
    for line in fp:
        if flag:
            parts = line.split(' ')
            parts = [ele for ele in parts if ele != '']
            parts = map(float, parts)
            score += parts
            if len(score) == 20:
                print "{",
                buff = ''
                for i in xrange(20):
                    if i == 0:
                        buff += "'{}': {:7.3f}".format(aa[i], score[i])
                    else:
                        buff += ",'{}': {:7.3f}".format(aa[i], score[i])
                buff += "},"
                print buff
                score = []
                flag = False
        if re.match(r'^I', line):
            flag = True
