import gzip
import re
from tensorflow.keras.models import model_from_json

from .. import Interval, Sequence, Variant, Orientation, Coordinates
from . import reverseComplement

from enum import Enum
class FileType(Enum):
    VCF = 1
    TSV = 2
    BED = 3

class IntervalIO():
    @staticmethod
    def getIntervals(file, format=FileType.BED):
        intervals = []
        with gzip.open(file, 'rt') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                lineSplit = line.split("\t")
                id = None
                if (len(lineSplit) >=4):
                    id = lineSplit[3]
                if (len(lineSplit) >=6 and lineSplit[5].strip() == "-"):
                    intervals.append(Interval(lineSplit[0], int(lineSplit[2]), int(lineSplit[1])+1, id))
                else:
                    intervals.append(Interval(lineSplit[0], int(lineSplit[1])+1, int(lineSplit[2]),id))
        return(intervals)



class VariantIO():
    @staticmethod
    def loadVariants(file, removeRefVariants=False, fileType=FileType.VCF):
        variants = []
        with gzip.open(file, 'rt') as f:
            for line in f:
                split_line = line.strip().split("\t")
                if line.startswith("#"):
                    continue
                if (fileType == FileType.VCF):
                    variant = Variant(split_line[0],int(split_line[1]),split_line[3],split_line[4])
                if (fileType == FileType.TSV):
                    variant = Variant(split_line[0],int(split_line[1]),split_line[2],split_line[3])
                if (not removeRefVariants or not variant.isRef()):
                    variants.append(variant)

        return(variants)

class SequenceIO():
    @staticmethod
    def readFasta(fastaFile, skipN=True, coordinates=Coordinates.ONEBASED_INCLUSIVE_INCLUSIVE):
        sequences = []
        if (fastaFile[-3:] == ".gz"):
            f = gzip.open(fastaFile, 'rt')
        else:
            f=  open(fastaFile, 'rt')

        id = None
        interval = None
        sequence = ""
        p = re.compile('^>(chr)?((1[0-9]|2[0-2]|X|Y|M|[0-9]):[0-9]+-[0-9]+){1}(\(([+-.])\))?(\s|_)?(.*)')
        for line in f:
            if line.startswith(">"):
                if id is not None:
                    sequences.append(Sequence(interval.contig, interval.start(), interval.end(), id, sequence))
                if ":" not in line: # WORKAROUND
                    if (line.count("_") > 4): # Skip unplaced contig sequences. This I have to do before! (e.g. chr16_KI270853v1_alt_525051_525350_neg_489)
                        id = None
                        continue
                    line = line.replace("_",":",1).replace("_",'-',1)
                match = p.search(line)

                interval = Interval.fromString("chr"+match.group(2), coordinates=coordinates)

                if (match.group(5) != None):
                    interval.orientation = Orientation(match.group(5))

                id = match.group(7)

                sequence = ""
            else:
                sequence += line.strip()
                if (skipN):
                    if ("N" in line):
                        id = None

        if id is not None:
            sequences.append(Sequence(interval.contig, interval.start(), interval.end(), id, sequence))
        f.close()
        return(sequences)

    @staticmethod
    def writeFasta(fastaFile,sequences):
        if (fastaFile[-3:] == ".gz"):
            f = gzip.open(fastaFile, 'wt')
        else:
            f=  open(fastaFile, 'wt')
        for sequence in sequences:
            f.write(sequence.getFastaString())
            f.write("\n")
        f.close()
    @staticmethod
    def getSequences(intervals, reference):
        sequences = []
        for interval in intervals:
            sequences.append(SequenceIO.readSequence(reference,interval))
        return(sequences)

    # Read a sequence from a pyfaidx Fasta file
    # 1 based and inclusive
    # returns sequence in upper case
    @staticmethod
    def readSequence(reference,interval, id=None):
        start = interval.start()
        end = interval.end()

        if (interval.isReverse()):
            start = interval.end()
            end = interval.start()

        sequence = reference[interval.contig][start-1:end].seq.upper()

        if (interval.isReverse()):
            sequence = reverseComplement(sequence)
        return(Sequence(interval.contig, interval.start(), interval.end(),  id, sequence))

class ModelIO():

    @staticmethod
    def loadModel(model,weights=None):
        json_file = open(model, 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        if (weights is not None):
            loaded_model.load_weights(weights)
        return(loaded_model)

    @staticmethod
    def saveModel(model,modelpath,weights=None):
        model_json = model.to_json()
        with open(modelpath, "w") as json_file:
            json_file.write(model_json)
        if (weights is not None):
            model.save(weights)
