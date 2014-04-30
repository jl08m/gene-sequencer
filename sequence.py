__author__ = 'Jia'

import itertools


class Sequence:

    #if more bases used, add to list
    BASES = ['A', 'C', 'G', 'T']
    NUM_BASES = len(BASES)

    def __init__(self, original_sequence=None):
        if original_sequence is None:
            raise StandardError('sequence cannot be none')

        self.sequence = original_sequence
        self.length = len(original_sequence)

    def generate_sequences(self, bases=1, spacer_length=0, spacer_position=1, start_position=1, end_position=None):
        """ Main method of class. Returns a list of permutated sequences based on the arguments.

        Generalized steps of method:
            The section that requires mutations will be the main-slice.
            The original sequence will be cut into pre-slice, main-slice, and post-slice.
            Generate a list of all permutations for the main-slice.
            Re-attach the pre- and post- slices to each of the permutated main-slices.
            return list of re-attached slices.

        Arguments:
        bases -- the number of bases that the main-slice will have.
        spacer_length -- the number of bases that a 'spacer' will have inside the main-slice.
        spacer_position -- where in the main-slice the spacer will start.
        start_position -- where in the original sequence the main-slice starts.
        end_position -- where in the original sequence the main-slice ends.

        """
        #set end to length of sequence
        if end_position is None:
            end_position = self.length

        permutated_sequences = set()

        #Sequence position starts at 1-n. Must fix to match proper array notation.
        start_position -= 1

        #Initialize slices.
        slice_start = start_position
        slice_end = slice_start + bases
        spacer_end = spacer_position + spacer_length

        while slice_end <= end_position:
            #get sequence parts before and after slice
            preslice_sequence = self.sequence[:slice_start]
            postslice_sequence = self.sequence[slice_end:]
            permutated_slices = []

            slice_sequence = self.sequence[slice_start:slice_end]

            #part of slice that gets permutated
            permutation_piece = slice_sequence[:spacer_position] + slice_sequence[spacer_end:]

            spacer = slice_sequence[spacer_position:spacer_end]

            #get permutations
            #FIXME: this part can be optimized. Duplicate permutations possible.
            temp = itertools.product(self.BASES, repeat=len(permutation_piece))
            permutations = set()
            for tmp in temp:
                permutations.add(str(tmp))

            #deprecated permutation part
            ##perms = permuteAll(permutation_piece)
            ##get rid of duplicate starting sequence in permutated list
            ##del perms[0]

            #reattach permutations with spacer
            for permutation in permutations:
                permutated_slices.append(permutation[:spacer_position] + spacer + permutation[spacer_position:])

            #put slice back into sequence parts
            for permutated_slice in permutated_slices:
                permutated_sequences.add(preslice_sequence + permutated_slice + postslice_sequence)

            slice_start += 1
            slice_end = slice_start + bases

        return permutated_sequences

    def _validate(self, bases, spacer_length, spacer_position, start, end):

        if isinstance(bases, int) or \
                isinstance(spacer_length, int) or isinstance(spacer_position, int) or \
                isinstance(start, int) or isinstance(end, int):
            raise TypeError("An argument is not an integer")

        if bases < 1:
            raise ValueError("Bases to change cannot be less than 1")
        if bases > self.length:
            raise ValueError("Bases to change cannot be greater than length of sequence")

        if spacer_length > (bases-2):
            raise ValueError("Spacer length exceeds allowed bases-2")
        if spacer_position + spacer_length >= bases:
            raise ValueError("Spacer position is invalid with given spacer lenght")

        if start < 1:
            raise ValueError("Start Position cannot be less than 1")
        if start > self.length:
            raise ValueError("Start Position cannot be greater than length of Sequence")

        if end < start:
            raise ValueError("End Position cannot be less than start position")
        if end > self.length:
            raise ValueError("End Position cannot be more than length of Sequence")

        return True

    def sortlist(sequence, sequence_list):
        """Returns a sequence list with the index sequence first"""
        index = sequence_list.index(sequence)
        retlist = sequence_list[index:] + sequence_list[:index]
        return retlist

