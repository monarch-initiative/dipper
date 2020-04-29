/*
View "mgd.mrk_summary_view"
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
    mt.name AS subtype,
    (((m.symbol || ', '::text) || m.name) || ', Chr '::text) || m.chromosome AS description,
    m.symbol AS short_description
   FROM acc_accession a
    join acc_accession a2 on a._object_key = a2._object_key
    join mrk_marker m on a._object_key = m._marker_key
    join mrk_types mt on m._marker_type_key = mt._marker_type_key 
  WHERE m._organism_key = 1 
    AND a._mgitype_key = 2 
    AND a.private = 0 
    AND a2._logicaldb_key = 1 
    AND a2._mgitype_key = 2 
    AND a2.prefixpart = 'MGI:'::text 
    AND a2.preferred = 1
  ;
