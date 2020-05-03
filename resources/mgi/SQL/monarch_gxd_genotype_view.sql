/*
First derived from View "mgd.gxd_genotype_view"

~ 5 seconds
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
   FROM gxd_genotype g
    join prb_strain s on g._strain_key = s._strain_key
    join acc_accession a on g._genotype_key = a._object_key 
    join voc_term vt on g._existsas_key = vt._term_key
    join acc_mgitype at on a._mgitype_key =  at._mgitype_key
  WHERE at.name = 'Genotype'
    AND a._logicaldb_key = 1 
    AND a.prefixpart = 'MGI:'::text 
    AND a.preferred = 1 
 ;
