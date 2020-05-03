/*
First derived from View "mgd.gxd_allelepair_view"

 ~ 5 seconds

*/
 SELECT a._allelepair_key,
    a._genotype_key,
    a._allele_key_1,
    a._allele_key_2,
    a._marker_key,
    a._mutantcellline_key_1,
    a._mutantcellline_key_2,
    a._pairstate_key,
    a._compound_key,
    a.sequencenum,
    m.symbol,
    m.chromosome,
    a1.symbol AS allele1,
    a2.symbol AS allele2,
    t1.term AS allelestate,
    t2.term AS compound
   FROM gxd_allelepair a
     JOIN mrk_marker m ON a._marker_key = m._marker_key
     JOIN all_allele a1 ON a._allele_key_1 = a1._allele_key
     JOIN all_allele a2 ON a._allele_key_2 = a2._allele_key
     JOIN voc_term t1 ON a._pairstate_key = t1._term_key
     JOIN voc_term t2 ON a._compound_key = t2._term_key
UNION
 SELECT a._allelepair_key,
    a._genotype_key,
    a._allele_key_1,
    a._allele_key_2,
    a._marker_key,
    a._mutantcellline_key_1,
    a._mutantcellline_key_2,
    a._pairstate_key,
    a._compound_key,
    a.sequencenum,
    m.symbol,
    m.chromosome,
    a1.symbol AS allele1,
    NULL::text AS allele2,
    t1.term AS allelestate,
    t2.term AS compound
   FROM gxd_allelepair a
     JOIN mrk_marker m ON a._marker_key = m._marker_key
     JOIN all_allele a1 ON a._allele_key_1 = a1._allele_key
     JOIN voc_term t1 ON a._pairstate_key = t1._term_key
     JOIN voc_term t2 ON a._compound_key = t2._term_key
  WHERE a._allele_key_2 IS NULL
;
