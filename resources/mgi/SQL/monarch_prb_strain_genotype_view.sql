/*
View "mgd.prb_strain_genotype_view"
*/
 SELECT s._straingenotype_key,
    s._strain_key,
    s._genotype_key,
    s._qualifier_key,
    ss.strain,
    gs.strain AS description,
    a.accid AS mgiid,
    t.term AS qualifier,
   FROM prb_strain_genotype s
    join prb_strain ss on s._strain_key = ss._strain_key 
    join acc_accession a on s._genotype_key = a._object_key 
    join gxd_genotype g on s._genotype_key = g._genotype_key 
    join prb_strain gs on g._strain_key = gs._strain_key 
    join voc_term t on s._qualifier_key = t._term_key 
  WHERE a._mgitype_key = 12 
    AND a._logicaldb_key = 1 

  ;
