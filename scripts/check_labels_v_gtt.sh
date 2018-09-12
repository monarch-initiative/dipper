#! /bin/bash

# check that lables found in the dipper code
# are aslo keys in the global translation table.
#
# surfaces mis spelled or missing mappings

IFS=':'
GTT=translationtable/GLOBAL_TERMS.yaml
TMPLABEL=/tmp/lables_found

egrep -Roh "self.globaltt\['[^']+'\]" dipper |\
	cut -f2 -d "'" | sort -u | sed 's|\(.*\)|^"\1": |g' > "${TMPLABEL}"

for lab in $(<"${TMPLABEL}"); do
	# echo "${lab}"
	if(! grep -qo "${lab}" "${GTT}"); then
			echo "${lab} not found in ${GTT}"
	fi
done

#rm "${TMPLABEL}"
