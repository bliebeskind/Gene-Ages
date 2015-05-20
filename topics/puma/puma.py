'''
Probabilistic Unigram Mixture Allocation - a variation on LDA where each document is drawn from but one topic.
This simpler clustering method may be more suitable to documents which are brief and sparse.
Keegan Hines
05/19/15
'''

import numpy.random as r
from collections import Counter
from itertools import groupby
from functools import reduce  

class model():
    def __init__(self):
        pass  
    
    def simData(self,D):
        '''
        Simulate some fake data for sanity checks. Data drawn from 3 topics, D documents, with vocabulary of size 4.
        Superbly inflexible method for synthetic data generation, but all we need for now.
        '''
        thetas = [ [.9,.1,0,0],[.25, .25, .25, .25],[0,.1,.1,.8]]
        docs = []
        for i in range(D):
            pr = thetas[r.choice([0,1,2])]
            c = Counter(r.choice(range(len(pr)), size= r.poisson(15),p = pr  ))
            docs.append( [c.get(s,0) for s in [0,1,2,3] ] )
        return docs
            
        
    def probability_doc_given_topics(self,doc,topic_probs):
        '''
        Given the parameters of a single topic, return likelihood of the doc given the parameters. 
        Input
            doc: List[Int]
            topics_probs: List[Double]
        '''
       
        return reduce(lambda a,b:a*b, [b**a for (a,b) in zip(doc,topic_probs)])
    
    def topicSample(self, doc, params):
        '''
        Resample the topic label of the doc.
        Inputs
            doc: List[Int]
            params
        '''
        probabilities = []
        for T in params:
            
            probabilities.append( self.probability_doc_given_topics(doc[1],T[1]) )
        probabilities = [p/float(sum(probabilities)) for p in probabilities]
        draw = r.choice(range(len(params)),1,p=probabilities)
        return (draw[0], doc[1])
    
    def parameterSample(self, allData):
        '''
        Given the document labels, resample the parameters for each topic.
        '''
        
        allData.sort(key = lambda a: a[0])
        grouped = [(k, [thing[1] for thing in g]) for k,g in groupby(allData, lambda s: s[0])]
        group_counts = [(group[0], reduce(lambda a,b: [a[i] + b[i] for i in range(len(a))]   ,group[1]) ) for group in grouped]
        group_probs = [(k, r.dirichlet(c).tolist()) for k,c in group_counts]
        return group_probs
        
        
    def run(self,corpus,num_topics, iterations):
        
        # initial data structure as list of (label, [counts])
        # ex [ (1, [0,1,0,8]),  (0, [2,5,0,0]), (1, [3,4,9,1])  ]
        
        labeledCorpus = map( lambda c: (r.choice(range(num_topics)) ,c) ,corpus)
        
        for i in range(iterations):
        
            # sample the parameters given the labels
            
            topicParameters = self.parameterSample(labeledCorpus)
            
            # sample the labels given the parameters
            labeledCorpus = map(lambda doc: self.topicSample(doc,topicParameters) ,labeledCorpus)
        return topicParameters
        
        
if __name__ == '__main__':
    p = model()
    corp = p.simData(1800)
    res = p.run(corp,num_topics=3,iterations=50)
    print res    