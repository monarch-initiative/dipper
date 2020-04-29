/*
View "mgd.gxd_genotype_view"
*/
 SELECT g._genotype_key,
    g._strain_key,
    g.isconditional,
    g.note,
    g._existsas_key,
    s.strain,
    a.accid AS mgiid,
    (('['::text || a.accid) || '] '::text) || s.strain AS displayit,
    vt.term AS existsas
   FROM gxd_genotype g,
    prb_strain s,
    acc_accession a,
    voc_term vt,
  WHERE g._strain_key = s._strain_key 
    AND g._genotype_key = a._object_key 
    AND a._mgitype_key = 12 
    AND a._logicaldb_key = 1 
    AND a.prefixpart = 'MGI:'::text 
    AND a.preferred = 1 
    AND g._existsas_key = vt._term_key
 ;
