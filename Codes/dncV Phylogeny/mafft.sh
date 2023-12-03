#!bin/bash

#required:
## Input: fasta file containing the amino acid sequences of dncVs

mafft  --retree 1 --clustalout --reorder dncV.fasta > dncv.nwk

echo done
