

def modify_sequences_by_variants(sequences, variants):
    for sequence in sequences:
        for variant in variants:
            sequence.replace(variant)
    return(sequences)


def reverseComplement(sequence):
    return(complement(reverse(sequence)))

def reverse(sequence):
    return(sequence[::-1])
def complement(sequence):
    out=''
    for i in sequence:
        if (i == 'A' or i == 'a'):
            out += 'T'
        elif (i =='T' or i == 't'):
            out += 'A'
        elif (i =='C' or i == 'c'):
            out += 'G'
        elif (i =='G' or i == 'g'):
            out += 'C'
        else:
            out += i
    return(out)
