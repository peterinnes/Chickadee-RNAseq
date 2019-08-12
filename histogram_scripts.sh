####getting a data table to make histogram of fixdif freqs in N-mt genes####
grep -wFf mito_genes_MitoCarta moch_vs_bcch/snpEff_genes_moch_bcch_BCCH_fixdif.txt | sort -k 8 -rn | tail -n +14 | cut -f 1,8 | sort -u > moch_vs_bcch/N-mt_genes_plt-hist

##getting data table for number of fixed differences in non-mito genes##
grep -wvFf mito_genes_MitoCarta moch_vs_bcch/snpEff_genes_moch_bcch_BCCH_fixdif_SORTED.txt | cut -f 1,5 > non-mt_genes_plt-hist

##add catagory column for mt, N-mt, and non-mt genes##
tail -n +2 non-mt_genes_plt-his | while read line; do sed 's/$/\tnon-mt/'; done > temp #problem: while loop is skipping the first line



##awk one liner to sum across rows, i.e. sum fix diffences from each impact category##
awk '{ for (i=2; i<=NF; i++) j+=$i; print j; j=0}' temp | paste temp - | less
