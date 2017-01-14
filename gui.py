import easygui as q

import os

r = q.diropenbox()
L = os.listdir(r)
choices_dict = {}
for item in L:
    choices_dict[item[:4]] =item

PROJECT = q.choicebox(msg='Cancer Project',choices=choices_dict.keys())

print PROJECT
PICKLE_INFILE = os.path.join(r, choices_dict[PROJECT])
print PICKLE_INFILE
