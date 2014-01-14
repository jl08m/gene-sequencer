

===============================================================================

GENERATE SEQUENCE SCRIPT: generate_sequences.py

It generates a list of sequences with all possible combinations given 
the parameters with no repeats.
Output stored in JSON format
	"sequence length" -- length of the sequence 
	"range start" -- start of the range, from the fasta file 
	"number of sequences" -- number of sequences in the list
    "gene id" -- the gene's id from the fasta (dna id)
	"chrom" -- chromosome from the fasta file
	"range end" -- the end of the range from the fasta file
	"sequence" -- actual sequence
	"key" -- identifying key number
	"mutations start" -- where the first mutation is
	"mutations" -- how many mutations total
	"mutations end" -- where the last mutation is

===============================================================================

PARAMETER RULES:

K-mers:
	the number of bases you wish to change.
		Example: Changing 2 bases means
		AAAAA->CAAAA->CCAAA->...->AAATT
	It must encompass the spacer length.
	
Spacer Length:
	How many bases you want to skip over.
		Example: 4kmer, 2 spacers means
		A**AA->C**TA->...->T**TA
	Where ** is the spacer of any base and it does not change.
		
Spacer Position:
	Where in the kmer grouping you want the spacer placed
		Example: 4kmer, 1 spacer, after position 2
		12345  12345       12345
		AA*AA->CA*GA->...->TT*TA
	
Starting and Ending position:
	Where in the whole input sequence you want the mutations to be generated
	
Most likely errors are spacer length and position being larger than the kmer 
or out of bounds.

===============================================================================

CONFIGURATION FILE : config.cfg

This is the file that all the scripts pull data from. There are many options
but do not change the option's name. The name is hardcoded so if you change it
the code most also be changed. 

===============================================================================

BATCH SCRIPT : batch.py

This script is the one that takes all the generated sequences and gives us data.
It generates the svm for each of the sequences and compares them against the
original sequence. This gives the correlation r-value that is used to sort
the final output list.

Output stored in JSON format sorted by descending r-value
	
	"sequence length" -- length of the sequence 
	"range start" -- start of the range, from the fasta file 
	"number of sequences" -- number of sequences in the list
    "gene id" -- the gene's id from the fasta (dna id)
	"chrom" -- chromosome from the fasta file
	"range end" -- the end of the range from the fasta file
	"base sequence" -- original sequence and its values
	"list":
		"sequence" -- actual sequence
		"key" -- identifying key number
		"mutations start" -- where the first mutation is
		"mutations" -- how many mutations total
		"mutations end" -- where the last mutation is
		"svmvals" -- an array of svm values of the sequence
		"rval" -- array of [r-value, p-value]

===============================================================================

SELECT SCRIPT : select_sequences.py

Generates a list of all the sequences that match the percentage cutoff point.
So a top and bottom of 5% would give a list of all sequences that fell within
the top 5% and bottom 5% of the full list.

	"sequence length" -- length of the sequence 
	"range start" -- start of the range, from the fasta file 
	"number of sequences" -- number of sequences in the list
    "gene id" -- the gene's id from the fasta (dna id)
	"chrom" -- chromosome from the fasta file
	"base sequence" -- original sequence and its values
	"list":
		"sequence" -- actual sequence
		"key" -- identifying key number
		"mutations start" -- where the first mutation is
		"mutations" -- how many mutations total
		"mutations end" -- where the last mutation is
		"svmvals" -- an array of svm values of the sequence
		"rval" -- array of [r-value, p-value]

===============================================================================

Create GFF : creategff.py

Generates the GFF files necessary. Only produces those GFF's in the selected
data file. Folder outputs are hardcoded so it is best if those were not renamed.

===============================================================================

R script : graphs.r

Run the R script after the creategff.py script is run to get a pdf of all the
plots comparing each sequence to the base sequence.

===============================================================================

CLEAN UP

All the GFF files must be moved/deleted after each run or else they will
accumulate and be used in the R script that creates the graphs.
Other files generated such as the list of sequences and selected data files
will be overwritten during the next run if the same filenames are used.

===============================================================================
