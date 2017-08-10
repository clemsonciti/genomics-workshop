#!/usr/bin/env python

import sys
import string
from collections import Counter
import argparse
import time
import multiprocessing

def merge_dicts(d1, d2):
    return dict(Counter(d1) + Counter(d2))

def word_count(fname):
    with open(fname) as f:
        data = f.read()

    data_without_punctuation = data.translate(None, string.punctuation)
    data_without_newlines = ' '.join(data_without_punctuation.split('\n'))
    data_without_newlines = data_without_newlines.lower()
    data_as_string_list = data_without_newlines.split(' ')

    c = Counter(data_as_string_list)
    time.sleep(5)
    return dict(c)

parser=argparse.ArgumentParser(description=
        """ Print word frequency for one or more text files.
        Removes all punctuation and converts words to lowercase
        before counting.
        """)
parser.add_argument('-j', type=int, metavar='N', default=1, help='Number of tasks to use')
parser.add_argument('files', type=str, help='Input files', nargs='+')

args = parser.parse_args()
files = args.files
jobs = args.j

results = []

pool = multiprocessing.Pool(jobs)
results = pool.map(word_count, files)
pool.close()
pool.join()

word_freqs = {}
for result in results:
    word_freqs = merge_dicts(word_freqs, result)

for word in word_freqs:
    print("{} {}".format(word_freqs[word], word))
