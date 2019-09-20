#!/usr/bin/env python2

from p4 import *

import argparse

parser = argparse.ArgumentParser(description='combining Supports')
parser.add_argument('-l', '--list', help='taxa list(alignemts in nexus format)', required=True)
parser.add_argument('-m', '--master', help='master tree(in nexus format)', required=True)
parser.add_argument('-s', '--second', help='second tree(in nexus format)', required=True)
parser.add_argument('-o', '--outfile', help='outfile', default='combinedTree.nex')
parser.add_argument('-r', '--root', help='root', default='manual')
args = parser.parse_args()
master = args.master
second = args.second
taxalist = args.list
out = args.outfile
root = args.root
# a.taxNames is the list we want.
var.doCheckForDuplicateSequences=False
read(taxalist)

a = var.alignments[0]

with open(master, 'r') as file :
  filedata = file.read()


filedata = filedata.replace('[&label=', '')
filedata = filedata.replace(']', '')
filedata = filedata.replace('[&R', '')  
# Write the file out again
with open(master, 'w') as file:
  file.write(filedata)




read(master)

var.punctuation = var.phylip_punctuation
tMB = var.trees[0]  # name the first

with open(second, 'r') as file :
  filedata = file.read()


filedata = filedata.replace('[&label=', '')
filedata = filedata.replace(']', '')
filedata = filedata.replace('[&R', '')  
# Write the file out again
with open(second, 'w') as file:
  file.write(filedata)


read(second)
tPAUP = var.trees[1] # name it

tMB.taxNames = a.taxNames
tPAUP.taxNames = a.taxNames

# Split keys are numerical versions of the 'dot-star' split notation.
# The same split on the two trees would have the same split key.
tMB.makeSplitKeys()
tPAUP.makeSplitKeys()

# Make a dictionary, so that we can fish out nodes in the paup tree
# given a split key.  Split keys are found on node branches, here
# n.br.
myDict = {}
for n in tPAUP.iterInternalsNoRoot():
    myDict[n.br.splitKey] = n

for nMB in tMB.iterInternalsNoRoot():
    # Given a split key in the mrbayes tree, we can find the
    # corresponding node in the paup tree, using the split key with
    # the dictionary.
    nPAUP = myDict.get(nMB.br.splitKey)
    # If there was none, then nPAUP is None
    if nPAUP:
        nMB.name = '%s/%s' % (nMB.name, nPAUP.name)
    else:
        nMB.name = '%s/-' % nMB.name
# rooting the tree   
if root == 'manual':
    pass
else:
    n = tMB.node(root).parent
    tMB.reRoot(n)

#print nMB.name
tMB.writeNewick(out)

