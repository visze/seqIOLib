import numpy as np


class Encoder():

    #this is set up for 1d convolutions where examples
    #have dimensions (len, num_channels)
    #the channel axis is the axis for one-hot encoding.
    @classmethod
    def one_hot_encode_along_channel_axis(cls, sequence):
        to_return = np.zeros((len(sequence),4), dtype=np.int8)
        cls._seq_to_one_hot_fill_in_array(to_return, sequence, 1)
        return to_return

    @classmethod
    def _seq_to_one_hot_fill_in_array(cls, zeros_array, sequence, one_hot_axis):
        assert one_hot_axis==0 or one_hot_axis==1
        if (one_hot_axis==0):
            assert zeros_array.shape[1] == len(sequence)
        elif (one_hot_axis==1):
            assert zeros_array.shape[0] == len(sequence)
        #will mutate zeros_array
        for (i,char) in enumerate(sequence):
            if (char=="A" or char=="a"):
                char_idx = [0]
            elif (char=="C" or char=="c"):
                char_idx = [1]
            elif (char=="G" or char=="g"):
                char_idx = [2]
            elif (char=="T" or char=="t"):
                char_idx = [3]
            elif (char=="Y" or char=="y"):
                char_idx = [1,3]
            elif (char=="R" or char=="r"):
                char_idx = [0,2]
            elif (char=="S" or char=="s"):
                char_idx = [1,2]
            elif (char=="W" or char=="w"):
                char_idx = [0,3]
            elif (char=="K" or char=="k"):
                char_idx = [2,3]
            elif (char=="M" or char=="m"):
                char_idx = [0,1]
            elif (char=="B" or char=="b"):
                char_idx = [1,2,3]
            elif (char=="D" or char=="d"):
                char_idx = [0,2,3]
            elif (char=="H" or char=="h"):
                char_idx = [0,1,3]
            elif (char=="V" or char=="v"):
                char_idx = [0,1,2]
            elif (char=="N" or char=="n"):
                continue #leave that pos as all 0's
            else:
                raise SequenceCharacterError("Unsupported character", char)
            if (one_hot_axis==0):
                zeros_array[char_idx,i] = 1
            elif (one_hot_axis==1):
                zeros_array[i,char_idx] = 1

class SequenceCharacterError(Exception):
    def __init__(self, message, character):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self,message)

        self.character = character

    def __str__(self):
        return("%s: %s" % (self.args[0], str(self.character)))
