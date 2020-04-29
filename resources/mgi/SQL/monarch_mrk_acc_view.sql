/*
View "mgd.mrk_acc_view"
*/
 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.private,
    a.preferred

    l.name AS logicaldb,
    m._organism_key
   FROM acc_accession a,
    acc_logicaldb l,
    mrk_marker m
  WHERE a._mgitype_key = 2 
    AND a._logicaldb_key = l._logicaldb_key 
    AND a._object_key = m._marker_key
 ;
