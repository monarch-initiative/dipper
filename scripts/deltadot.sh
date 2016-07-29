#! /bin/bash

# deltadot.sh  one.dot two.dot > delta.dot

# Use to illuminate the commonality and differences between dot files created 
# by 'ntriples2dot.awk' but might work for other dot files 

# makeup names for temporary files
edge1="$(mktemp -q)" || exit 1
edge2="$(mktemp -q)" || exit 2

# isolate and order edges, here including edge label
# but strip the edge counts included by ntriples.awk as they are non structural
grep "\->" "$1"| sed 's/ ([[:digit:]]*)"];/"];/g' | sort > "${edge1}"
grep "\->" "$2"| sed 's/ ([[:digit:]]*)"];/"];/g' | sort > "${edge2}"

# Could color nodes as well but their state is easily determined by their edges

# combine graph names in their given order
label1="$(grep '^label=' "${1}" | cut -f2 -d \")"
label2="$(grep '^label=' "${2}" | cut -f2 -d \")"

# output
echo "digraph {"
echo "rankdir=LR;"
echo 'charset="utf-8";'
echo "LITERAL [shape=record];"
echo 'labelloc="t";'
echo "label=\"${label1} <=> ${label2}\""
### 	partition the edges 
# green are in both dot files
comm -12 "${edge1}" "${edge2}"|sed 's/"];/", color="green"];/g'
# red are not in second dot
comm -23 "${edge1}" "${edge2}"|sed 's/"];/", color="red"];/g'
# blue are not in first dot
comm -13 "${edge1}" "${edge2}"|sed 's/"];/", color="blue"];/g'
echo  -e "\n}\n"

# cleanup (typically /tmp/tmp.WHATEVER)
rm "${edge1}" "${edge2}"
