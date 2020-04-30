/*
View "mgd.prb_strain_acc_view":
*/
 SELECT a._accession_key,
    a.accid,
    'ComplexReigon:' || a._object_key as mgi_internal,
    al.name AS logical_name
   FROM acc_accession a
    join acc_logicaldb al on a._logicaldb_key = al._logicaldb_key
    join mrk_types mt on a._mgitype_key = mt._mgitype_key
  WHERE  mt.name = 'Complex/Cluster/Region' 
    AND a.private = 0
    AND a.preferred = 1
  ;
