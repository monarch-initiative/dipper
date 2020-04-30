/*
View "mgd.mrk_acc_view"
*/
 SELECT a._accession_key,
    a.accid,
    at.name || ':' || a._object_key as mgi_internal,
    l.name AS logicaldb,
    m._organism_key
   FROM acc_accession a
    join acc_logicaldb l on a._logicaldb_key = l._logicaldb_key 
    join mrk_marker m on a._object_key = m._marker_key
    join acc_mgitype at on a._mgitype_key =  at._mgitype_key
  WHERE at.name = 'Marker'
  	AND a.private = 0
  	and a.preferred = 1
 ;
