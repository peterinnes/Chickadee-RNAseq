#!/bin/bash

set -exuo pipefail

if [ $# -lt 1 ]
then
    echo "
    Searches for candidate genes in the snpEff_genes output, which has counts of variants per gene of varying impact level. This script will sum total number of variants per gene.
    Input is snpEff_genes output, and a list of gene symbols.
    
    [-s] Path to snpEff_genes summary file
    [-g] Path to gene list. One symbol per line.

    "

else
    while getopts s:g: option
    do
    case "${option}"
    in


