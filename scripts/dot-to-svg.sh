#!/usr/bin/env bash

for dotfile in `find . -type f -name "*.dot"`
do
    fname=`basename $dotfile .dot`
	dot -Tsvg $dotfile > docs/img/${fname}.svg
done
