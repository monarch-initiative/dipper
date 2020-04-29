/*
View "mgd.prb_strain_acc_view":
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
    l.name AS logicaldb
   FROM acc_accession a,
    acc_logicaldb l
  WHERE a._mgitype_key = 10 
  AND a._logicaldb_key = l._logicaldb_key
  ;
