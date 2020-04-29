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
    a.private,
    a.preferred,
    a._createdby_key,
    a._modifiedby_key,
    a.creation_date,
    a.modification_date,
    a2.accid AS mgiid,
    t.term AS subtype,
    (al.symbol || ', '::text) || al.name AS description,
    al.symbol AS short_description
   FROM acc_accession a,
    acc_accession a2,
    all_allele al,
    voc_term t
  WHERE a._mgitype_key = 11 
    AND a.private = 0 
    AND a._object_key = a2._object_key 
    AND a2._logicaldb_key = 1 
    AND a2._mgitype_key = 11 
    AND a2.prefixpart = 'MGI:'::text 
    AND a2.preferred = 1 
    AND a._object_key = al._allele_key 
    AND al._allele_type_key = t._term_key
;

