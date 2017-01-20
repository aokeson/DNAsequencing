import sys
import getopt
import random

def main():

    seqlength, readlength, errorperc = parsecommand()

    Basic(seqlength, readlength)
    Repeat(seqlength, readlength)
    Errors(seqlength, readlength, errorperc)
    Coverage(seqlength, readlength)

def parsecommand():
    #Define variables
    seqlength = 0
    readlength = 0
    errorperc = 0

    #get command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 's:r:e:')
    except getopt.GetoptError as err:
        print(err)
        print('Command line arguments were incorrect. Please rerun with correct arguments.')
        sys.exit(1)

    #put the arguments into the right variables
    for opt, arg in opts:
        if opt == '-s':
            seqlength = arg
        if opt == '-r':
            readlength = arg
        if opt == '-e':
            errorperc = arg

    # Check that mandatory arguments were entered
    if seqlength == 0 or readlength == 0:
        print("Arguments not entered. Please run again with correct arguments: -f, -g, -s.")
        sys.exit(1)


    return seqlength, readlength, errorperc


#Generates simple sequential reads with overlap to check if program works
def Basic(seqlength, readlength):
    sequence = ''
    counter = 0
    overlap = 10  #Edit this to change overlap between reads
    readlen = int(readlength)
    seqlen  = int(seqlength)
    fastafile = 'Basic.fasta'
    nucl = ['A','C','G','T']
    outstring =''

    #Generate sequence
    for i in range(seqlen):
        sequence += random.choice(nucl)
    print('>Basic')
    print(sequence)

    fp = open(fastafile,"w")  #open file

    for j in range(0,seqlen-readlen,readlen-overlap):
        outstring = '>read_' + str(counter) + '\n'
        fp.write(outstring)
        outstring = sequence[j:j+readlen] + '\n'
        fp.write(outstring)
        counter += 1

    outstring = '>read_' + str(counter) + '\n'
    fp.write(outstring)
    outstring = sequence[seqlen-readlen:]
    fp.write(outstring)

    return

#Generates reads with repeated regions
def Repeat(seqlength, readlength):
    sequence = ''
    counter = 0
    overlap = 80  #Edit this to change overlap between reads
    readlen = int(readlength)
    seqlen  = int(seqlength)
    fastafile = 'Repeat.fasta'
    nucl = ['A','C','G','T']
    outstring =''

    #Generate sequence
    for i in range(int(seqlen/3)):
        sequence += random.choice(nucl)

    sequence += sequence[-100:] #repeat last 100 nucleotides
    sequence += sequence[-150:] #repeat last 150 nucleotides

    for i in range(seqlen-len(sequence)):
        sequence += random.choice(nucl)
    print('>Repeat')
    print(sequence)

    fp = open(fastafile,"w")  #open file

    for j in range(0,seqlen-readlen,readlen-overlap):
        outstring = '>read_' + str(counter) + '\n'
        fp.write(outstring)
        outstring = sequence[j:j+readlen] + '\n'
        fp.write(outstring)
        counter += 1

    outstring = '>read_' + str(counter) + '\n'
    fp.write(outstring)
    outstring = sequence[seqlen-readlen:]
    fp.write(outstring)

    return

#Generates reads with errors
def Errors(seqlength, readlength, errorperc):

    sequence = ''
    counter = 0
    errcounter = 0
    overlap = 10  #Edit this to change overlap between reads
    readlen = int(readlength)
    seqlen  = int(seqlength)
    fastafile = 'Error.fasta'
    nucl = ['A','C','G','T']
    outstring =''

    #Generate sequence
    for i in range(seqlen):
        sequence += random.choice(nucl)
    print('>Errors')
    print(sequence)

    fp = open(fastafile,"w")  #open file

    for k in range(10):   #Coverage amount
        for j in range(0,seqlen-readlen,readlen-overlap):
            outstring = '>read_' + str(counter) + '\n'
            fp.write(outstring)
            seqlist = list(sequence[j:j+readlen])
            for a in range(len(seqlist)):
                if(float(errorperc) > random.random()):
                    seqlist[a] = random.choice(nucl)
                    errcounter += 1

            outstring = ''.join(seqlist) + '\n'
            fp.write(outstring)
            counter += 1

        outstring = '>read_' + str(counter) + '\n'
        fp.write(outstring)
        seqlist = list(sequence[seqlen-readlen:])
        for a in range(len(seqlist)):
            if(float(errorperc) > random.random()):
                seqlist[a] = random.choice(nucl)
                errcounter += 1
        outstring = ''.join(seqlist) + '\n'
        fp.write(outstring)
        counter += 1

    print('Error Counter: ',errcounter)
    return

#Generates random reads from the sequence so some areas will have more reads than others
def Coverage(seqlength, readlength):
    sequence = ''
    counter = 0
    coverage = 10           #change this to alter coverage
    readlen = int(readlength)
    seqlen  = int(seqlength)
    fastafile = 'Coverage.fasta'
    nucl = ['A','C','G','T']
    outstring =''

    #Generate sequence
    for i in range(seqlen):
        sequence += random.choice(nucl)
    print('>Coverage')
    print(sequence)

    fp = open(fastafile,"w")  #open file
    for k in range(coverage*int(seqlen/readlen)):
        j = random.randint(0,seqlen)
        outstring = '>read_' + str(counter) + '\n'
        fp.write(outstring)
        if(j+readlen<seqlen):
            outstring = sequence[j:j+readlen] + '\n'
            fp.write(outstring)
        else:
            outstring = sequence[j-readlen:j] + '\n'
            fp.write(outstring)
        counter += 1

    return

main()