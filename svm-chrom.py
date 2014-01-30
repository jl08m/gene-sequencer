__author__ = 'Jia'
#!/bin/env python
# FILE: svm-chrom.py
# AUTHOR: William Stafford Noble and Shobhit Gupta
# CREATE DATE: 2008-03-24
import sys
import string
import logging
import filecmp
usage = """USAGE: svm-chrom.py <chromosome> <model>\n"""


#############################################################################
def __make_kmer_list(k, alphabet):
    # Base case.
    if (k == 1):
        return (alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return ([])

    # Error case.
    if (k < 1):
        logging.error("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in __make_kmer_list(k - 1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])

    return (return_value)

##############################################################################
def __make_upto_kmer_list(k_values, alphabet):
    # Compute the k-mer for each value of k.
    return_value = []
    for k in k_values:
        return_value.extend(__make_kmer_list(k, alphabet))

    return (return_value)

##############################################################################
def __find_revcomp(sequence, revcomp_dictionary):
    # Save time by storing reverse complements in a hash.
    if (revcomp_dictionary.has_key(sequence)):
        return (revcomp_dictionary[sequence])

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if (letter == "A"):
            return_value = return_value + "T"
        elif (letter == "C"):
            return_value = return_value + "G"
        elif (letter == "G"):
            return_value = return_value + "C"
        elif (letter == "T"):
            return_value = return_value + "A"
        elif (letter == "N"):
            return_value = return_value + "N"
        else:
            sys.stderr.write("Unknown DNA character (%s)\n" % letter)
            sys.exit(1)

    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value

    return (return_value)

##############################################################################
def __frequency_normalize(k_values, vector, kmer_list):
    # Initialize all vector lengths to zeroes.
    vector_lengths = {}
    for k in k_values:
        vector_lengths[k] = 0

    # Compute sum or sum-of-squares separately for each k.
    num_kmers = len(kmer_list)
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        vector_lengths[kmer_length] += count

    # Divide through by each sum.
    return_value = []
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        vector_length = vector_lengths[kmer_length]
        if (vector_length == 0):
            return_value.append(0)
        else:
            return_value.append(float(count) / float(vector_length))

    return (vector_lengths, return_value)

##############################################################################
def __update_sequence_vector(seq_add,
                             seq_rem,
                             vector,
                             discriminant,
                             hyperplane,
                             revcomp_dictionary,
                             kmer_indices,
                             sum_vector_lengths):
    ## Adjust for the Added kmers
    for i in range(1, len(seq_add) + 1):
        kmer = seq_add[len(seq_rem) - i:]
        rev_kmer = __find_revcomp(kmer, revcomp_dictionary)

        # Take the kmer with A or C
        if (cmp(kmer, rev_kmer) > 0):
            kmer = rev_kmer

        # Initial determinant value
        d_init = vector[kmer_indices[kmer]] * hyperplane[kmer_indices[kmer]]

        # Change in vector
        vector[kmer_indices[kmer]] = ((vector[kmer_indices[kmer]] \
                                       * sum_vector_lengths[len(kmer)]) + 1) \
                                     / sum_vector_lengths[len(kmer)]

        # Final determinant value
        d_after = vector[kmer_indices[kmer]] * hyperplane[kmer_indices[kmer]]

        # Adjust the discriminant
        discriminant = discriminant - d_init + d_after

    ## Adjust for the removed kmers
    for i in range(1, len(seq_rem) + 1):
        kmer = seq_rem[0:i]
        rev_kmer = __find_revcomp(kmer, revcomp_dictionary)

        # Take the kmer with A or C
        if (cmp(kmer, rev_kmer) > 0):
            kmer = rev_kmer

        # Initial determinant value
        d_init = vector[kmer_indices[kmer]] * hyperplane[kmer_indices[kmer]]

        # Change in vector
        vector[kmer_indices[kmer]] = ((vector[kmer_indices[kmer]] \
                                       * sum_vector_lengths[len(kmer)]) - 1) \
                                     / sum_vector_lengths[len(kmer)]

        # Final determinant value
        d_after = vector[kmer_indices[kmer]] * hyperplane[kmer_indices[kmer]]

        # Adjust the discriminant
        discriminant = discriminant - d_init + d_after

    return (vector, discriminant)

##############################################################################
def __make_sequence_vector(sequence,
                           revcomp_dictionary,
                           k_values,
                           alphabet,
                           kmer_list):
    num_ns = 0
    # Make an empty counts vector.
    kmer_counts = {}

    # Iterate along the sequence.
    for k in k_values:
        seq_length = len(sequence) - k + 1
        for i_seq in range(0, seq_length):

            # Extract this k-mer.
            kmer = sequence[i_seq: i_seq + k]

            # Store the count in the version that starts with A or C.
            rev_kmer = __find_revcomp(kmer, revcomp_dictionary)
            if (cmp(kmer, rev_kmer) > 0):
                kmer = rev_kmer

            # Increment the count.
            if (kmer_counts.has_key(kmer)):
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1

    # Build the sequence vector.
    sequence_vector = []
    for kmer in kmer_list:
        if (kmer_counts.has_key(kmer)):
            sequence_vector.append(kmer_counts[kmer])
        else:
            sequence_vector.append(0)

    # Normalize it
    (sum_vector_lengths, return_value) = __frequency_normalize(k_values,
                                                               sequence_vector,
                                                               kmer_list)

    return (sequence.count('N'), sum_vector_lengths, return_value)

##############################################################################
# MAIN
##############################################################################
f = open('temp/test_a375/svm.txt', 'w')

# Parse the command line.
if (len(sys.argv) != 3):
    sys.stderr.write(usage)
    sys.exit(1)
chrom_filename = sys.argv[1]
model_filename = sys.argv[2]

# Open the chromosome file.
if (chrom_filename == "-"):
    chrom_file = sys.stdin
else:
    chrom_file = open(chrom_filename, "r")

# Read the first line.
chrom_id = chrom_file.readline()
chrom_id = chrom_id[1:-1]

# Read the chromosome into memory.
chromosome = ""
for line in chrom_file:
    line = line[:-1].upper()
    chromosome = chromosome + line
sys.stderr.write("Read %d characters in %s.\n" % (len(chromosome), chrom_id))

# Read the model.
model_file = open(model_filename, "r")

# Bias is on the penultimate line.
prev_line = ""
line_num = 1
for line in model_file:

    # C value is on the fifth line.
    if (line_num == 5):
        regularizer = float(line.split()[0])
        sys.stderr.write("regularizer=%g\n" % regularizer)

    bias = prev_line
    prev_line = line
    line_num += 1
bias = float(bias.split()[0])
sys.stderr.write("bias=%g\n" % bias)

# Hyperplane coordinates are on the last line.
hyperplane = line.split()
for index in range(0, len(hyperplane)):
    hyperplane[index] = float(hyperplane[index])

# The alphabet.
alphabet = "ACGT"

# Make a list of all values of k.
k_values = range(1, 7)

# Make a list of all k-mers.
kmer_list = __make_upto_kmer_list(k_values, alphabet)
sys.stderr.write("Considering %d kmers.\n" % len(kmer_list))

# Set up a dictionary to cache reverse complements.
revcomp_dictionary = {}

# Use lexicographically first version of {kmer, revcomp(kmer)}.
new_kmer_list = []
for kmer in kmer_list:
    rev_kmer = __find_revcomp(kmer, revcomp_dictionary)
    if (cmp(kmer, rev_kmer) <= 0):
        new_kmer_list.append(kmer)
kmer_list = new_kmer_list;
sys.stderr.write("Reduced to %d kmers.\n" % len(kmer_list))

# Store kmer indices
kmer_indices = {}
index = 0
while index < len(kmer_list):
    kmer_indices[kmer_list[index]] = index
    index += 1

# Verify that the model width is the same.
if (len(hyperplane) != len(kmer_list)):
    sys.stderr.write("Error: %d k-mers and %d hyperplane coordinates.\n"
                     % (len(kmer_list), len(hyperplane)))
    sys.exit(1)

# Extract the first subsequence.
position = 0
sequence = chromosome[position:position + 50]

# Compute the sequence vector.
(num_ns, sum_vector_lengths, vector) = \
    __make_sequence_vector(sequence,
                           revcomp_dictionary,
                           k_values,
                           alphabet,
                           kmer_list)

# Compute the initial discriminant.
discriminant = bias
for index in range(0, len(vector)):
    discriminant += vector[index] * hyperplane[index]
f.write("%g\n" % (discriminant / regularizer))

## Special case for speedup: Do you know the score for a sequence of 50Ns?
dis_50n = -245324
if num_ns == 50:
    dis_50n = discriminant

# Boolean: Do we have a computed vector?
computed_vector = True

# Extract remaining overlapping subsequences.
for position in xrange(1, len(chromosome) - 49):
    # The sequence
    sequence = chromosome[position:position + 50]

    # Extract new leading and trailing k-mers.
    seq_add = chromosome[position + 44:position + 50]
    seq_rem = chromosome[position - 1:position + 5]

    # Check if the added position is N or not
    unchanged = True
    if chromosome[position - 1] == 'N':
        num_ns -= 1
        unchanged = False
    if chromosome[position + 49] == 'N':
        num_ns += 1
        unchanged = False

    #print sequence, num_ns, seq_add, seq_rem, "hi", chromosome[position-1], chromosome[position+49]

    if ((num_ns == 0) and (unchanged)):
        # Update the vector and the corresponding discriminant.
        (vector, discriminant) = __update_sequence_vector \
                (seq_add,
                 seq_rem,
                 vector,
                 discriminant,
                 hyperplane,
                 revcomp_dictionary,
                 kmer_indices,
                 sum_vector_lengths)
    else:
        if (num_ns == 50):
            if dis_50n == -245324:
                (num_ns, sum_vector_lengths, vector) = \
                    __make_sequence_vector(sequence,
                                           revcomp_dictionary,
                                           k_values,
                                           alphabet,
                                           kmer_list)
                # Compute the initial discriminant.
                discriminant = bias
                for index in range(0, len(vector)):
                    discriminant += vector[index] * hyperplane[index]

                # Assign the score for 50Ns future use
                dis_50n = discriminant
            else:
                discriminant = dis_50n
        else:
            (num_ns, sum_vector_lengths, vector) = \
                __make_sequence_vector(sequence,
                                       revcomp_dictionary,
                                       k_values,
                                       alphabet,
                                       kmer_list)
            # Compute the initial discriminant.
            discriminant = bias
            for index in range(0, len(vector)):
                discriminant += vector[index] * hyperplane[index]

    # Tell the user what's happening.
    if (position > 0) and (position % 1000 == 0):
        sys.stderr.write("Position=%d.\n" % position)

    f.write("%g\n" % (discriminant / regularizer))
f.close()

# Local variables:
# mode: python
# py-indent-offset: 2
# End:
