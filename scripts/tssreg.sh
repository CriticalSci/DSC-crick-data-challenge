#!/bin/bash
## Configuration modified from dgoguerra https://gist.github.com/dgoguerra/9206418

### File name
readonly PROGNAME=$(basename $0)
### File name, without the extension
readonly PROGBASENAME=${PROGNAME%.*}
### File directory
readonly PROGDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
### Arguments
readonly ARGS="$@"
### Arguments number
readonly ARGNUM="$#"

usage() {
    echo "Script description"
    echo
    echo "Usage: $PROGNAME -i <file> -p <file> -g <file> -f <file> -l <integer> -r <integer> [options]..."
    echo
    echo "Options:"
    echo
    echo "  -h, --help"
    echo "      This help text."
    echo
    echo "  -i <genes.bed>, --infile <genes.bed>"
    echo "      BED6 File with gene coordinates. If \"-\", stdin will be used instead."
    echo
    echo "  -p <peaks.bed>, --infile <peaks.bed>"
    echo "      BED6 File with peak coordinates. If \"-\", stdin will be used instead."
    echo
    echo "  -g <reference.genome>, --genome <reference.genome>"
    echo "      Genome Index File with chromosome sizes. If \"-\", stdin will be used instead."
    echo
    echo "  -f <reference.fasta>, --fasta <reference.fasta>"
    echo "      Genome FASTA File. If \"-\", stdin will be used instead."
    echo
    echo "  -l <integer>, --left <integer>"
    echo "      5' Offset from TSS (positive or negative). If \"-\", stdin will be used instead."
    echo
    echo "  -r <integer>, --right <integer>"
    echo "      3' Offset from TSS (positive or negative). If \"-\", stdin will be used instead."
    echo
    echo "  --"
    echo "      Do not interpret any more arguments as options."
    echo
}

while [ "$#" -gt 0 ]
do
    case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
    -i|--infile)
        infile="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    -p|--peaks)
        peaks="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    -g|--genome)
        genome="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    -f|--fasta)
        fasta="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    -l|--left)
        left="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    -r|--right)
        right="$2"

        # Jump over <file>, in case "-" is a valid input file 
        # (keyword to standard input). Jumping here prevents reaching
        # "-*)" case when parsing <file>
        shift
        ;;
    --)
        break
        ;;
    -*)
        echo "Invalid option '$1'. Use --help to see the valid options" >&2
        exit 1
        ;;
    # an option argument, continue
    *)  ;;
    esac
    shift
done

## Commands
prefix="${infile%.*}_${left}${right}"

awk -v offset5=$left -v offset3=$right '
BEGIN {
    OFS = "\t";
    #offset5 = -2000;
    #offset3 = 0;
    }
NR == FNR {
    chr[$1] = $2;
    next
}
{
    if ($1 in chr) {
        lim = chr[$1];
        if ($6 == "+") {
            left = $2 + offset5;
            right = $2 + offset3
        } else {
            left = $3 - offset3;
            right = $3 - offset5
        };
        if (left < 0) {
            left = 0
        } else if (right > lim) {
            right = lim
        } else {
                $2 = left;
                $3 = right
            };
        print
    } else {
        print "error "$1" not in genome"
    }
}' $genome $infile | \
bedtools intersect -a - -b $peaks | \
sort -k1,1 -k2,2n | \
bedtools getfasta -name -s -fi $fasta -bed - | \
awk '/^>/{$0=$0"_"(++i)}1' >&1 #> "${prefix}.fasta"


