# Gene-Ages

## A consensus approach to estimating gene ages for model organisms

This repository contains scripts, data, and ipython notebooks in support of our manuscript:

**"Towards Consensus Gene Ages"** 
Benjamin J. Liebeskind, Claire D. McWhite, and Edward M. Marcotte

### How old is my gene?

If you study a certain gene or gene family, you might be interested in knowing how many other organisms share orthologs of that gene. Or maybe you want to annotate genomic datasets with gene ages to get an idea of how deep in evolutionary time different pathways were assembled. If so, look no further!

We took orthology calls from 13 popular orthology inference algorithms and estimated consensus gene ages for a variety of model eukaryotes. Because orthology inference is notoriously difficult, we also annotated our datasets with various error terms so that you can propagate uncertainty through your downstream analyses.

### Organisms with gene-age information

You can find consensus tables for the following organisms in the Main/ directory. They are named main_\<UniprotID\>.csv

| **Common Name** | **Uniprot ID** | 
| --------------- | -------------- |
| Anopheles gambiae (Mosquito) | ANOGA | 
| Bos taurus (Cattle) | BOVIN | 
| Branchiostoma floridae (Lancelet) | BRAFL | 
| Caenorhabditis elegans (Worm) | CAEEL | 
| Candida albicans | CANAL | 
| Canis lupus familiaris (Dog) | CANFA | 
| Gallus gallus (Chicken) | CHICK | 
| Ciona intestinalis (Tunicate) | CIOIN | 
| Cryptococcus neoformans | CRYNJ | 
| Danio rerio (Zebrafish) | DANRE | 
| Drosophila melanogaster (Fly) | DROME | 
| Homo sapiens (Human) | HUMAN | 
| Ixodes scapularis (Tick) | IXOSC | 
| Macaca mulatta (Rhesus macaque) | MACMU | 
| Monosiga brevicollis (Choanoflagellate) | MONBE | 
| Monodelphis domestica (Opossum) | MONDO | 
| Mus musculus (Mouse) | MOUSE | 
| Nematostella vectensis (Sea anemone) | NEMVE | 
| Neurospora crassa | NEUCR | 
| Ornithorhynchus anatinus (Platypus) | ORNAN | 
| Pan troglodytes (Chimp) | PANTR | 
| Phaeosphaeria nodorum | PHANO | 
| Rattus rattus (Rat) | RAT | 
| Saccaromyces cerevisiae (Budding yeast) | YEAST | 
| Schistosoma mansoni (Blood fluke) | SCHMA | 
| Schizosaccharomyces pombe (Fission yeast) | SCHPO | 
| Sclerotinia sclerotiorum | SCLS1 | 
| Takifugu rubripes (Pufferfish) | TAKRU | 
| Ustilago maydis (Corn smut/Huitlacoche | USTMA | 
| Xenopus tropicalis (Frog) | XENTR | 
| Yarrowia lipolytica | YARLI |

### Age-categories

These files contain the following information:
- First, they contain a distribution over gene ages estimated by 13 orthology algorithms
after trimming errors (see manuscript for details). This distribution can be viewed as
a posterior.
- Next, they contain the mode age of the gene
- The next 6 columns contain details on various sources of error in the age-call

#### Error statistics

| Name | Description |
| ---- | ----------- | 
| NumDBsContributing | How many databases/algorithms contribute to final estimate. More is better |	
| NumDBsFiltered | How many databases/algorithms were trimmed out. Less is better |
| entropy | Shannon's entropy over final age-call distribution. Lower is better |
| NodeError	| Average patristic distance between age calls before filtering. Lower is better |
| Bimodality | How bimodal the age call is (see manuscript). Lower is better | 
| HGT_flag | Whether or not this gene was flagged as being a recent horizontal gene transfer |

### Replicating the analysis

We provide scripts and ipython notebooks if you're interested in replicating the analysis or
running again with some different parameters. You should run the code using scripts in 
CannedScripts/. Here's a flowchart to show how the scripts in CannedScripts/, some of the
ipyton notebooks in Notebooks/, and the output files in Data/ are all related:

![Flowchart](https://github.com/bliebeskind/Gene-Ages/blob/master/pics/FlowChart.png)


And because you can now add emojis to GitHub markdown, I will
:black_joker: