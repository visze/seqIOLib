from pyfaidx import Fasta
import math, gzip
import copy

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

from enum import Enum
class Orientation(Enum):
    FORWARD = '+'
    REVERSE = '-'
    UNKNOWN = '.'

class Coordinates(Enum):
    ZEROBASED_INCLUSIVE_EXCLUSIVE = 1
    ZEROBASED_INCLUSIVE_INCLUSIVE = 2
    ONEBASED_INCLUSIVE_EXCLUSIVE = 3
    ONEBASED_INCLUSIVE_INCLUSIVE = 4

class VariantType(Enum):
    SNV = 1
    INSERTION = 2
    DELETION = 3

class Position:
    def __init__(self,contig,position):
        self.contig = contig
        self.position = position
    def __str__(self):
        return("%s:%d" % (self.contig,self.position))

class Variant(Position):
    def __init__(self,contig,position, ref, alt):
        Position.__init__(self,contig,position)
        self.ref = ref
        self.alt = alt
        if len(ref) == len(alt):
            self.type = VariantType.SNV
        elif (len(ref) > len(alt)):
            self.type = VariantType.DELETION
        elif (len(ref) < len(alt)):
            self.type = VariantType.INSERTION

    def __str__(self):
        return("%s:%d%s>%s"% (self.contig,self.position,self.ref,self.alt))

    # checks if a variant is just the reference sequence. Meaning ref and als is the same or alt is defined as a dot
    def isRef(self):
        return(self.ref == self.alt or self.alt == ".")

# 1 based
class Interval(Position):
    def __init__(self,contig,start,end,id=None):
        self.id = id
        if (end < start):
            self.orientation = Orientation.REVERSE
            Position.__init__(self,contig,end)
            self.length = start-end+1
        else:
            self.orientation = Orientation.FORWARD
            Position.__init__(self,contig,start)
            self.length = end-start+1

    @classmethod
    def fromString(cls, string, coordinates=Coordinates.ONEBASED_INCLUSIVE_INCLUSIVE):
        region = string.replace("-",":").split(":")
        contig = region[0]
        start = int(region[1])
        end = int(region[2])
        if (coordinates == Coordinates.ZEROBASED_INCLUSIVE_EXCLUSIVE):
            if (start > end):
                end += 1
            else:
                start += 1
        elif (coordinates == Coordinates.ZEROBASED_INCLUSIVE_INCLUSIVE):
            end += 1
            start += 1
        elif (coordinates == Coordinates.ONEBASED_INCLUSIVE_EXCLUSIVE):
            if (start > end):
                start += 1
            else:
                end += 1
        return cls(contig, start, end)

    def __str__(self):
        return("%s:%d-%d" % (self.contig,self.start(),self.end()))

    def start(self):
        if (self.isReverse()):
            return(self.position+self.length-1)
        else:
            return(self.position)

    def end(self):
        if (self.isReverse()):
            return(self.position)
        else:
            return(self.position+self.length-1)

    def isReverse(self):
        return(self.orientation == Orientation.REVERSE)

    def contains(self,position):
        return(self.contig == position.contig and self.position <= position.position and self.position+self.length > position.position)

    def include(self,interval):
        return(self.contains(Position(interval.contig,interval.start())) and self.contains(Position(interval.contig,interval.end())))

    def tiling(self, length=300, shift=1):

        if (length < shift):
            exit("Shift must have same or smaller size then length")

        if ((shift * (self.length//shift+1) + (length % shift)) % self.length == 0):
            missing = 0
        elif (length % shift != 0):
            exit("length have to be dividable without reset by shift")
        elif (self.length > length):
            missing=self.length % shift
        elif (self.length < length):
            missing=length - self.length
        else:
            exit("Interval length with shifts did not fit n times into interval")

        new_start = self.position - math.ceil(missing/2)
        new_end = self.position + self.length-1 + math.floor(missing/2)
        tiling_start = new_start
        tiling_end = tiling_start + length - 1

        intervals = []
        while tiling_end <= new_end:
            if (self.isReverse()):
                intervals.append(Interval(self.contig, tiling_end, tiling_start))
            else:
                intervals.append(Interval(self.contig, tiling_start, tiling_end))

            tiling_start += shift
            tiling_end = tiling_start + length -1
        return(intervals)

class Sequence(Interval):
    def __init__(self,contig,start,end,id,sequence):
        Interval.__init__(self,contig,start,end,id)

        if (self.orientation == Orientation.REVERSE):
            self.sequence=reverseComplement(sequence)
        else:
            self.sequence=sequence

    def __str__(self):
        return(self.getFastaString()[0:80])

    def getSubSequence(self,interval):
        if (self.include(interval)):

            startPos=interval.position-self.position
            newSequence = self.sequence[startPos:startPos+interval.length]
            result = Sequence(self.contig,interval.position, interval.position + interval.length-1,
                        self.id, newSequence)
            if((self.isReverse() and interval.isReverse) or (not self.isReverse() and interval.isReverse())):
                result.orientation = Orientation.REVERSE
            return(result)
        return None

    def replace(self, variant):
        if (self.contains(variant)):
            ref_len = len(variant.ref)
            # alt_len = len(variant.alt)
            start = variant.position - self.position
            self.sequence = self.sequence[0:start]+variant.alt+self.sequence[start+ref_len:]
        # if (start == 0):
        #     self.position = self.position +
        # else:
        #     self.length = self.length + alt_len-ref_len

    def getSequence(self, orientation=None):
        if (orientation == None):
            orientation=self.orientation
        if (orientation == Orientation.REVERSE):
            return(reverseComplement(self.sequence))
        else:
            return(self.sequence)

    def saturationMutagensis(self,start=1,end=0):
        start = start -1
        if (end <= 0):
            end = len(self.sequence) -1
        else:
            end = end - 1
        satMut = [copy.deepcopy(self)]
        variants = []
        for i, ref in enumerate(self.getSequence(Orientation.FORWARD)):
            if (i >= start and i<= end):
                for alt in ["A","C","G","T"]:
                    if (ref != alt):
                        variant = Variant(self.contig, self.position+i, ref, alt)
                        newSequence = copy.deepcopy(self)
                        newSequence.replace(variant)
                        satMut.append(newSequence)
                        variants.append(variant)
        return(satMut,variants)

    def getFastaString(self):
        start = self.start()
        end = self.end()
        if (self.isReverse()):
            end = self.start()
            start = self.end()

        id=""
        if (self.id is not None and len(self.id)>0):
            id = "_%s" % self.id
        return ">%s:%d-%d(%s)%s\n%s" % (self.contig,start,end,self.orientation.value,id,self.getSequence())



# reads a fasta sequence from the reference.
# 1) if the sequence is longer than the length it tiles it up into n sequences shifteted
#   by the shift parameter. Sequence will be extended in equal length left an right to fit the tilingOfRegion
# 2) if the sequence is smaller than the length: sequence is extended left and right to fit the length (centered)
# one based, both inclusive


def idForTiling(position, maxPosition, contig, start, end):
    id = "%s_tile%s/%d" % (intervalToString(contig, start, end), str(position).zfill(len(str(maxPosition))), maxPosition)
    return id


def tileRegionAndWriteFastaEntries(fout, reference, id, contig, start, end, length=300, shift=1, variants=None):
    tile_result = tilingOfRegion(reference, contig, start, end, length, shift)
    for i, tile_tuple in enumerate(tile_result):
        new_id = "%s_%s" % (id, idForTiling(i+1, len(tile_result), contig, start, end))
        sequence = Sequence(tile_tuple[0], tile_tuple[1], tile_tuple[2], new_id, tile_tuple[3])
        if (variants is not None):
            for variant in variants:
                sequence.replace(variant)
        writeFasta(fout, sequence)
