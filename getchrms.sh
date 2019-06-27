#! /bin/bash

fasta=$1

grep -E "^>" $1 | cut -d" " -f1 | tr -d '>'
