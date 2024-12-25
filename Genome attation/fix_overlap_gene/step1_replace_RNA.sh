sed -E 's/^([^[:space:]]+\s+[^[:space:]]+\s+)(C_gene_segment|V_gene_segment|primary_transcript|transcript)(\s+.*)$/\1mRNA\3/g;
        s/^([^[:space:]]+\s+[^[:space:]]+\s+)RNA(\s+.*)$/\1ncRNA\2/g' ../../chicken.v23.merge.longest_isoform.gff > step1.gff
