/*
View "mgd.all_summary_view"

*/

 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.preferred,
    a2.accid AS mgiid,
    t.term AS subtype,
    (al.symbol || ', '::text) || al.name AS description,
    al.symbol AS short_description
   FROM acc_accession a,
    join acc_accession a2 on a._object_key = a2._object_key 
    join all_allele al on a._object_key = al._allele_key 
    join voc_term t on al._allele_type_key = t._term_key
    join acc_mgitype at on a._mgitype_key  =  at._mgitype_key
    WHERE at.name = 'Allele' 
      AND a.private = 0 
      AND a._mgitype_key = a2._mgitype_key
      AND a2._logicaldb_key = 1 
      AND a2.prefixpart = 'MGI:'::text 
      AND a2.preferred = 1 
;

