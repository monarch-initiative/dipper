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

    l.name AS logicaldb,
    m._organism_key
   FROM acc_accession a
    join acc_logicaldb l on a._logicaldb_key = l._logicaldb_key 
    join mrk_marker m on a._object_key = m._marker_key
  WHERE a._mgitype_key = 2 
  	AND a.private = 0
  	and a.preferred = 1
 ;
