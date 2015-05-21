'''
Topic modeling with Ben's protein evolution data 
Keegan Hines
05/20/15
'''

import puma
from collections import namedtuple


ProteinComplex = namedtuple("ProteinComplex", ["Name", "Counts"])

def dataIngest(filename):
    '''
    Read in the data, parse it.
    Return a list of those ProteinComplex objects, and also the column header names for later use.
    '''
    with open(filename, "r") as inFile:
        txt = inFile.read()
    lines = txt.split("\n") 
    columnNames = lines[0].split("\t")[1:]
    dataset =  [ProteinComplex(line.split("\t")[0],[int(i) for i in line.split("\t")[1:]]) for line in lines[1:-1] ] # some rough python...
    return dataset, columnNames

if __name__ == "__main__":
    dataset, columnNames = dataIngest("../fixtures/corum_phylome_ages.txt")
    vocab_dictionary = dict(zip(range(len(columnNames)), columnNames ))
    
    corpus = map(lambda l: l[1],dataset)
    p = puma.model()
    result = p.run(corpus,num_topics=4, iterations= 20)
    
    print result