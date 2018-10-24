#!/bin/sh

# fps=15
examples_dir="examples/"
for d in $examples_dir* ; do
    # echo "$d"
    exmpl=$(echo $d| cut -d'/' -f 2)
    echo $exmpl
    python -u pgmt.py -t examples/$exmpl/$exmpl.tree -l examples/$exmpl/$exmpl-coords.txt -o output/$exmpl > output/$exmpl.log 2>&1
    sh ./gen_video.sh output/$exmpl
done
