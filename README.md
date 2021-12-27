## Welcome to Splinter (v1.0)

Splinter is a simple command line script intended to be used as a full coverage mutant library generator. It was originally created as a quick and nasty script for generating mutant DNA libraries. I could not find any applications on the internet that do this, so here's mine.

The research in [the lab where Iused to work] (http://yokobayashilab.net), focuses on probing the secondary structure of functional ribozymes, by assaying full coverage single and double replacement mutant libraries of some wildtype sequence of interest. See [Kobori, et al. 2015] (https://www.ncbi.nlm.nih.gov/pubmed/27461281) for an example of this research. This is a similar approach to that used in [Double Mutant Cycle Analysis] (http://www.sciencedirect.com/science/article/pii/S1359027896000569) for protein structure probing. 

## How Splinter Works 
Splinter runs in [Python3] (https://www.python.org/downloads/) from the command line. Just download the file mutagen.py and save it to your working directory. Then run it from the command line as follows:

```markdown
bioinfo:home$ python /home/path-to-working-directory/mutagen.py

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

etc...
```
Follow the command line prompts, and you're bound to come right.

### Dependencies
```markdown
[Numpy1.1] (http://www.numpy.org) or higher.
[pandas0.19] (http://pandas.pydata.org) or higher.
```
### Input

It takes from the user a single wildtype DNA sequence, with bases from non-mutable regions in lowercase and bases from the mutable regions in uppercase:

For example, your input/wildtype sequence is:

acACGTgt

All mutant bases in lowercase will be fixed, whereas the uppercase bases will be mutated. So for example, the input above would generate a single deletion library of four sequences for full coverage of the mutable region:

```markdown
ID	        Sequence	Mutation
Wildtype	acACGTgt	
0	        AC_CGTGT	A3, 
1	        ACAC_TGT	G5, 
2	        ACA_GTGT	C4, 
3	        ACACG_GT	T6, 
```
And similarly for the double deletion, triple deletion, single replacement and double replacement mutant libraries. All libraries are output to an excel file, which most manufacturers (e.g. [IDT] (http://sg.idtdna.com/site)) will accept to produce your oligonucleotides.

### Contact 

I have much to learn, so feel free to email me if you have any questions, comments, concerns and/or suggestions about the code. 

Crys

larvalanobium(at)gmail.com
