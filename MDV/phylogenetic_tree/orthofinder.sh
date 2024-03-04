#!/bin/bash
orthofinder -t 36 \
            -f /home/binz/MDV/result/OrthoFinder/proteins \
            -M msa \
            -S diamond \
            -a 24 \
            -A mafft \
            -T raxml-ng \

