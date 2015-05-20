'''
Topic modeling with Ben's protein evolution data 
Keegan Hines
05/19/15
'''

from collections import namedtuple
from gensim import models
import numpy as np

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
    dataset =  [ProteinComplex(line.split("\t")[0],[int(i) for i in line.split("\t")[1:]]) for line in lines[1:] ] # some rough python...
    return dataset, columnNames
    
def toSparse(l):
    '''
    Transform a dense list of counts into a sparse representation of (index, value) tuple list.
    example [0,1,6,0] becomes [(1,1), (2,6)]
    '''
    arr = np.array(l)
    return zip(arr.nonzero()[0], arr[arr>0])


if __name__ == "__main__":
    dataset, columnNames = dataIngest("./fixtures/corum_phylome_ages.txt")
    
    vocab_dictionary = dict(zip(range(len(columnNames)), columnNames ))
    K=4
    
    corpus = map(lambda s: toSparse( s.Counts ) ,dataset)
    
    lda= models.ldamodel.LdaModel(corpus, id2word=vocab_dictionary,num_topics=K)
    
    with open("lda_result","w") as writeFile:
        s = lda.print_topics(K, len(columnNames))
        [writeFile.write(topic + "\n") for topic in s]