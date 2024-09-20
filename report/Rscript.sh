#!/usr/bin/env bash


set -eo pipefail
module purge
module load mhipgx/R_Bioconductor/4.1.2_3.14


main() {
    Rscript "$@"
}

main "$@"
