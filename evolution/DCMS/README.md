# DCMS
```
/usr/bin/Rscript DCMS_linux.r -h
Usage: DCMS_202506.r [options]
This script performs positive and negative selection analysis using the DCMS method.


Options:
        --input-fst=CHARACTER
                Path to the input FST file (required)

        --input-pi1=CHARACTER
                Path to the input pi1 file for high group (required)

        --input-pi2=CHARACTER
                Path to the input pi2 file for low group (required)

        --species=CHARACTER
                Species to analyze: chicken, pig, cattle, sheep, dog [default: chicken]

        --out-prefix=CHARACTER
                Prefix for output files [default: output]

        --output-folder=CHARACTER
                Folder where output files will be saved [default: .]

        -h, --help
                Show this help message and exit

```
