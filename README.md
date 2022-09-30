
# Codon Usage Bias for Evolution (C.U.B.E)


###  Description
A codon is a set of three nucleotides in a sequence of DNA or RNA that corresponds to a specific amino acid. There are codon usage is degenerate as 64 different codons code for 20 amino acids. Thus multiple codons can code for same amino acids. All the codons are not used equally to express these amino acids. This preference to use a specific codon for an  amino acid is called Codon Usage Bias (CUB). Different evolutionary forces shape the CUB and they in turn affect the translation and gene expression levels making CUB very important. There are several statistical measures of CUB like Relative synonymous codon usage (RSCU), Neutrality plot analysis, Parity plot analysis, GC content, Effective number of codons etc. The Codon Usage Bias for Evolution (C.U.B.E) focusses on calculating the Relative synonymous codon usage (RSCU), GC content Neutrality plot analysis, Parity plot analysis of given sequences. It also gives the scatter plot for Neutrality plot analysis, Parity plot analysis and the extent of mutational pressure causing the CUB

### Workflow
The C.U.B.E reads a user-given input FASTA file consisting of coding sequences obtained from databases like GenBank. The code opens and reads the sequences and determines all the codons and their frequencies. This frequency is used to determine the RSCU values that are written in the output file names output.txt. The RSCU is the ratio of the observed frequency to the expected frequency for each synonymous codon. The code then calculates the total GC content and partial GC content of each sequence. The total GC content of each sequence is also written in the output file. The partial GC content is further used for Neutrality Analysis(GC12 vs GC3) and Parity Analysis (A3/AT3 vs G3/GC3). This data is plotted in two separate plots named “NeutralityGraph” and “ParityGraph”. The slope of the regression line in the Neutrality plot gives the extent of mutational pressure.

### Packages 
The required modules for statistical analysis and plotting the graphs are:
gonum.org/v1/gonum/stat
gonum.org/v1/plot”

### Usage
Run the file as 
./Project inputFileName.txt (for Mac)

### Results 
This program gives extent of mutational pressure for given coding sequences

### Applications
Translation dynamics
Gene expression levels
Codon Optimization




