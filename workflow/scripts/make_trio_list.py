#!/usr/bin/env python3

import glob
import fileinput
import os
import sys

index = []
mom = []
dad = []

for line in sys.argv[1]):
    if ',TE' in line and not 'Experiment Name' in line:
        for x in range(1, 5):
            if "Trio%s-Index" % x in line.split(',')[2]:
                index.append(line.split(',')[0])
            if "Trio%s-Foralder" % x in line.split(',')[2]:
                if "_M_" in line.split(',')[2]:
                    dad.append(line.split(',')[0])
                if "_K_" in line.split(',')[2]:
                    mom.append(line.split(',')[0])

original_stdout = sys.stdout # Save a reference to the original standard output

with open('sys.argv[2]', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    for i,j,k in zip(index,dad,mom):
        print(i,j,k)

sys.stdout = original_stdout # Reset the standard output to its original value
