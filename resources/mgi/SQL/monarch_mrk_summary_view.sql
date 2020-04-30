/*
View "mgd.mrk_summary_view"
*/
 SELECT a._accession_key,
    a.accid,
    al.name as logical_name,
    mt.name || ':' || a._object_key as mgi_internal,
    a.preferred,
    a2.accid AS mgiid,
    (((m.symbol || ', '::text) || m.name) || ', Chr '::text) || m.chromosome AS description,
    m.symbol AS short_description
   FROM acc_accession a
    join acc_accession a2 on a._object_key = a2._object_key
    join mrk_marker m on a._object_key = m._marker_key
    join mrk_types mt on m._marker_type_key = mt._marker_type_key 
    join mgi_orginism mo on m._organism_key = mo._organism_key
    join acc_logicaldb al on a._logicaldb_key = al._logicaldb_key
  WHERE mt.name = 'Marker'
  	AND mo.latinname = 'Mus musculus/domesticus' 
    AND a.private = 0 
    AND a2._logicaldb_key = 1 
    AND a2._mgitype_key = 2 
    AND a2.prefixpart = 'MGI:'::text 
    AND a2.preferred = 1
  ;
