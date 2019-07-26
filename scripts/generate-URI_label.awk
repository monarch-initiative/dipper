#! /usr/bin/awk -f
# ${DIPPER}/scripts/generate-URI_label.awk \
# ${DIPPER}/dipper/curie_map.yaml ${DIPPER}/translationtable/GLOBAL_TERMS.yaml

# Reconstitutes the label an URI has based on Dipper's mapping files.
# note: this relies on uninforced yaml quoting disipline

# in: curie_map.yaml    'pfx': 'buri'
NR == FNR && !/^#/{
	FS="'"
	a[$2]=$4
}

# in: GLOBAL_TERMS.yaml  "label": "curie",
NR != FNR && !/^#/{
	FS="\""
	split($4, t,":");
	if(t[1] in a)
		print a[t[1]] t[2] "\t" $2
}
# out: URI label
