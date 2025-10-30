#!/bin/bash

for fst in *.windowed.weir.fst; do
    if [[ ! -f "$fst" ]]; then
        continue
    fi
    
    prefix=$(basename "$fst" .windowed.weir.fst)
    
    pi1="${prefix}_high.windowed.pi"
    pi2="${prefix}_low.windowed.pi"
    
    if [[ -f "$pi1" && -f "$pi2" ]]; then
        echo "Processing: $fst with $pi1 and $pi2"
        
        /usr/bin/Rscript ~/project/01_evolution/script/DCMS_202506.r \
            --input-fst "$fst" \
            --input-pi1 "$pi1" \
            --input-pi2 "$pi2" \
            --species chicken \
            --out-prefix "$prefix" \
            --output-folder ./
    else
        echo "Skipping $fst: missing $pi1 or $pi2"
    fi
done

echo "All processing completed."
