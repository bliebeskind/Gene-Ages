from LECA.consensus import consensus_ages
import sys

### This program will create the consensus (mode) age calls 
### by trimming databases that oversplit co-orthologous groups.


############# User input #######################

INFILE = "binAges_human.csv"
LDORESULTS = "LDO_results.p"
FALSEPOSITIVES = "falsePos.p"
AGES = ['Cellular_organisms',
 'Euk+Bac',
 'Euk_Archaea',
 'Eukaryota',
 'Opisthokonta',
 'Eumetazoa',
 'Vertebrata',
 'Mammalia']

############ Don't change #######################

for line in consensus_ages(INFILE,AGES,LDORESULTS,FALSEPOSITIVES):
    print line
