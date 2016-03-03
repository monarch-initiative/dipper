#! /bin/bash

# Requires
# pip3 install pycallgraph
# [yum|apt-get] install graphviz



# relative path to where the output goes
outdir="./callgraph"

rm -fr ${outdir}

for src in `find . -type f -name "*.py" -exec echo \{\} \\;`; do 

	oldpth="${src%/*}"			    # drop all from last slash on
	basefn="${src##.*/}"			# drop all up to & inc last slash
	newpth="${outdir}${oldpth#.}"	# drop leading dot	
	rtname="${basefn%.py}"			# drop ".py" suffix

#    echo "src:    ${src}"

	if  [ ! -d ${newpth} ] ; then
		mkdir -p ${newpth}
#		echo "mkdir:  ${newpth}"
	fi
# 	echo "oldpth: ${oldpth}"
#	echo "basefn: ${basefn}"
#	echo "newpth: ${newpth}"
#	echo "rtname: ${rtname}"
#	echo ""


	echo "pycallgraph graphviz --output-file=\"${newpth}/${rtname}.png\" -- \"${src}\""
	echo ""
done
