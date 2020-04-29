/*
View "mgd.bib_acc_view"
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
    l.name AS logicaldb
  FROM acc_accession a,
    acc_logicaldb l
  WHERE a._mgitype_key = 1 
    AND a._logicaldb_key = l._logicaldb_key
;

