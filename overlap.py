######################################################################################
######################################################################################
#
#   File       : overlap.py
#
#   Purpose    : This program performs de novo assembly on next-generation
#			    sequencing data. The assembler uses a greedy overlap algorithm
#                to find the shortest common superstring (SCS) from the sequencing
#			    reads. Output is the set of SCS contigs pieced together by the 
#                assembly algorithm.
#
#   Developers : Alex Okeson
#			    Lee Korshoj
#                Karthik Handady
#                Andrew Gordon
#
#                CSCI 4314/5314, MCDB 5314 - Spring 2016
#                Final Project Application
#
######################################################################################
#
#   Functions:
#
#   main(argv) - parses, checks inputs, and calls functions determining SCS
#
#   usage() - helps user input arguments correctly
#
#   suffixPrefixMatch(seq1, seq2, overlap) - quantifies the overlap between 2 sequences
#
#   combine(matrix, seqNames, overlap) - performs the greedy SCS algorithm
#
#   readInput(inFile) - reads sequences and names into matrix
#
######################################################################################
# 
#   References: Lee Korshoj, Jan 2016
#               Korshoj_HW1.py
#               Formatting of the header content
#
#               David Knox (david.knox@colorado.edu), Jan 2016
#               Sample_Command_Parsing.py
#               Parsing command line arguments
#
#			   Ben Langmead
#			   Johns Hopkins University
#               Assembly & Shortest common superstring Lecture Slides
#               Algorithm idea and outline
#
######################################################################################

import sys, getopt, os.path

######################################################################################
#
#  Usage function 
#
#     Helps if the user is unfamiliar with the program or enters incorrect arguments.
#
#     Sample command: python overlap.py -f filename.FASTA -l 3
#
######################################################################################

def usage():
    
	# Print the usage of the application
    print 'python overlap.py -f <filename> -l <fraction>'
    print "-f specifies the FASTA file containing sequencing reads"
    print "-l specifies the smallest allowed fraction matching in the overlap regions (0.0 < l <= 1.0)"

######################################################################################
#  
#  Main application function
#
#    Parses and checks the input arguments.
#
#    Runs 'readInput' fuction to parse input file into a list of lists of the form
#    [[name, sequence], [name, sequence], ...] where name and sequence are strings.
#
#    Scores the prefix of each sequence compared to the suffix of every other sequence
#    using 'suffixPrefixMatch' and stores the value in a matrix.
#
#    Runs 'combine' function to merge the strings using the matrix of scores. The 
#    set of SCS contigs are returned.
#      
######################################################################################

def main(argv):
	
	inputFile = '' # The input file name
	overlap = 0 # The overlap length
	
	try:
		# Go through the input parameters to make sure there are -f and -l options
		opts, args = getopt.getopt(argv,"hf:l:")
	# If parsing goes wrong, print an error and exit
	except getopt.GetoptError:
		print "\nERROR! Command line parsing problem. Correct usage:"
		usage()
		sys.exit(2)
	for opt, arg in opts:
		# Process each input parameter
		if opt == '-h':
			print "\nCorrect usage is:"
			usage()
			sys.exit()
		if opt == '-f':
			if os.path.isfile(arg):
				inputFile = arg
			# If the input file given does not exist, print an error message and exit the program
			else:
				print "\nERROR! Input file must exist. Correct usage:"
				usage()
				sys.exit(2)
		if opt == '-l':
			overlap = float(arg)
	# If one of the arguments is missing, print an error message and exit the program
	if inputFile == '':
		print "\nERROR! An input file (-f) must be specified. Correct usage:"
		usage()
		sys.exit(2)
	if overlap <= 0 or overlap > 1.0:
		print "\nERROR! A smallest allowable fraction matching (-l) must be specified. Correct usage:"
		usage()
		sys.exit(2)

	# Print which FASTA file is being used
	print "\nUsing FASTA file " + inputFile

     # Run the 'readInput' function to format the input file
     # Store the returned list in variable sequences
	sequences = readInput(inputFile)

	# Initialize the overlap score matrix
	overlapMatrix = [[0 for x in range(len(sequences))] for x in range(len(sequences))]
	# Find the overlap score
	for i in range(len(sequences)):
		for j in range(len(sequences)):
			if i == j:
				overlapMatrix[i][j] = 0
			else:
				overlapMatrix[i][j] = suffixPrefixMatch(sequences[i][1], sequences[j][1], overlap)

	# Put all the sequence names in a list
	# Order of the list is the same as the order of the columns in the overlap matrix
	seqNames = ['' for x in range(len(sequences))]
	for i in range(len(sequences)):
		seqNames[i] = sequences[i][1]
	
	# Run the 'combine' algorithm to assemble the DNA subsequences into contigs
	output = combine(overlapMatrix, seqNames, overlap)
 
    # Output the contigs for the user
	for i in range(len(output)):
		print ">Result" + str(i+1) + " " + str(len(output[i]))
		print output[i]

######################################################################################
#
#  Prints a matrix by line for debugging and printing the original distance matrix
#
######################################################################################

def printMatrix(matrix):
	print 
	for m in matrix:
		print m

######################################################################################
#
#  Overlap determination function
#
#    Finds the length of the longest suffix of 'seq1' of length greater than or equal
#    to 'overlap' that matches a prefix of 'seq2'. If there is no overlap meeting this
#    criteria, the score is 0.
#
#    Inputs: The suffix sequence, 'seq1', the prefix sequence, 'seq2', and the length of the
#            minimum overlap length, 'overlap'.
#
#    Output: The overlap length of the prefix and suffix sequences.
#
######################################################################################

def suffixPrefixMatch(seq1, seq2, overlap):
    
    # Shorten the prefix sequence to the length of the suffix sequence if needed
    if len(seq1) > len(seq2):
        seq1=seq1[-len(seq2):]
    
    # Translate 'seq1' and 'seq2' across each other and count matches at each location
    for i in range(len(seq1)):
        matches=0
        for j in range(len(seq1)-i):
            if seq1[j+i] == seq2[j]:
                matches += 1
        # Stop when overlap meets the user-specified criteria
        if matches/float(len(seq1)-i) >= overlap and (len(seq1)-i) > 1:
            break         
    # Overlap score added to the overlap matrix is the number of nucleotides overlapping
    if len(seq1)-i == 1:
        return 0
    else:
        return len(seq1)-i 

######################################################################################
#
#  Greedy algorithm for finding SCS function
#
#    Finds the maximum of the overlap matrix and then combines those two sequences.
#    Then, the overlap matrix is recalculated based on the new sequence. The function
#    recurses until there is only 1 sequence left or there is no sequences that
#    overlap. All SCS contigs are returned.
#
#    Inputs: The overlap matrix, 'matrix', of the current sequences, column names,
#            'seqNames', of the sequence names in the same column order as matrix, and 
#            the minumum allowable overlap, 'overlap'.
#
#    Output: The assembled DNA contigs.
#
######################################################################################

def combine(matrix, seqNames, overlap):
    
	# If there is only 1 sequence left, this is the assembled sequence
	if len(seqNames) == 1:
		return seqNames
	else:
		# Initialize the maximum overlap matrix and the index of the maximum
		maxSoFar = 0
		index = [0,0]
		# Go through each element of the overlap matrix to find the maximum element
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				# If the element is greater than the current maximum, set it to be the new max
				if matrix[i][j]>maxSoFar:
					maxSoFar = matrix[i][j]
					index = [i,j]
		
		# If there is no overlapping sequence left, return the longest sequence that has been found
		if maxSoFar == 0:
			'''longestSeq = [0,0]
			for i in range(len(seqNames)):
				if len(seqNames[i]) > longestSeq[0]:
					longestSeq = [len(seqNames[i]),i]'''
			return seqNames

		# Find the combined sequence based on the highest overlap score
		newSeq = seqNames[index[0]][:-maxSoFar]+seqNames[index[1]]
		# New sequence is now where one of the old sequences was
		seqNames[index[1]] = newSeq
		# Remove the old sequence
		seqNames.remove(seqNames[index[0]])
		# If the sequence removed was in a column before where the new sequence was placed, update the index of the new sequence
		if index[1]>=index[0]:
			index[1] = index[1]-1

		# Remove the matrix row of the combined sequence
		matrix.remove(matrix[index[0]])
		# Remove the matrix column of the combined sequence
		for m in matrix:
			del m[index[0]]
		# Recalculate the overlap scores of the new sequence compared to every other sequence
		for i in range(len(matrix)):
			# Keep the diagonal zero
			if i != index[1]:
				# Recalculate the row
				matrix[index[1]][i] = suffixPrefixMatch(seqNames[index[1]], seqNames[i], overlap)
				# Recalculate the column
				matrix[i][index[1]] = suffixPrefixMatch(seqNames[i], seqNames[index[1]], overlap)

		# Recurse to combine two more sequences
		return combine(matrix,seqNames,overlap)

######################################################################################
#
#  Parsing input file function
#
#    Finds the beginning of each sequence and stores the sequence name in one
#    list and the corresponding sequence string in another list. The output combines
#    both lists.
#
######################################################################################

def readInput(inFile):
    
	f = open(inFile, 'r')
	sequenceNames = []
	sequences = []
	sequenceSoFar = []
	output = []
 
	for line in f.xreadlines():
		line = line.strip('\n')
		if line != '' and line[0] == '>':
			# Separate each sequence header by spaces and store the name (string before the first space)
			name = line[1:].split()
			sequenceNames.append(name[0])
			# Store the sequence string corresponding to the previous sequence in the file
			# Ensure the corresponding name and sequence stored with the same index
			if sequenceSoFar != []:
				sequences.append(''.join(sequenceSoFar))
				sequenceSoFar = []
		else:
			# Collect all the lines of a sequence in sequenceSoFar
			sequenceSoFar.append(line)
   
	# Add the last sequence string into sequences
	sequences.append(''.join(sequenceSoFar))
	# Make a list of lists with the names and sequences
	# Correct for capitalization inconsistencies
	for i in range(len(sequences)):
		output.append([sequenceNames[i],sequences[i].upper()])
  
	return output

######################################################################################
#
#  Calls the 'main' function to run the program
#
######################################################################################

if __name__ == "__main__":
	main(sys.argv[1:])
 
######################################################################################
######################################################################################