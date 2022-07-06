import keras
from keras.models import Sequential, Model, load_model
from keras import backend as K
import isolearn.keras as iso
import tensorflow as tf
import os
import pandas as pd
import numpy as np
import re
from aparent.predictor import *


##read in all sequences
##directly read in fasta file
def auto_open(fname,mode='rt',compression=None):
    """ Automatically create correct file object for
        normal, gzipped and bzipped files, so file objects are passed through
        the correct (de-)compressors.
    """
    ext=os.path.splitext(fname)[1].lower()
    if ext=='.gz' or compression=='gzip':
        if ext!='.gz':
            fname +='.gz'
        return gzip.open(fname,mode)
    else:
        return open(fname,mode)

def read_fasta_sequences(fname):
    with auto_open(fname) as fin:
        hd = None
        seq = None
        for line in fin:
            if line.startswith('>'):
                if hd is not None:
                    yield hd,"".join(seq)
                hd = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip().upper())
        if hd is not None:
            yield hd,"".join(seq)

fastafile="/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/denovo_hela_bed.TES205.fa"

model_path = '/binf-isilon/sandelin/people/mengjun/tools/aparent/saved_models/legacy_models/aparent_theano_legacy_30_31_34_padded.h5'

aparent_model = load_model(model_path)
aparent_encoder = get_aparent_encoder()
seqs = []
ids= []
for hd, seq in read_fasta_sequences(fastafile):
    seqs.append(seq)
    ids.append(hd)

##To do the prediction
print("predicting...")
iso_pred, cut_pred = aparent_model.predict(x=aparent_encoder(seqs))

logiso=logit(iso_pred)
res=pd.DataFrame(logiso, columns=["PAS_score"], index=ids)
print(res)
res.to_csv("/binf-isilon/sandelin/people/mengjun/Exosome_ML/data/PAS_score_denovo_TES205.csv", index=True, header=True)
