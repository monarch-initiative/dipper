/*
View "mgd.prb_strain_genotype_view"
*/
 SELECT s._straingenotype_key,
    s._strain_key,
    s._genotype_key,
    s._qualifier_key,
    s._createdby_key,
    s._modifiedby_key,
    s.creation_date,
    s.modification_date,
    ss.strain,
    gs.strain AS description,
    a.accid AS mgiid,
    t.term AS qualifier,
   FROM prb_strain_genotype s,
    prb_strain ss,
    acc_accession a,
    gxd_genotype g,
    prb_strain gs,
    voc_term t
  WHERE s._strain_key = ss._strain_key 
    AND s._qualifier_key = t._term_key 
    AND s._genotype_key = a._object_key 
    AND a._mgitype_key = 12 
    AND a._logicaldb_key = 1 
    AND s._genotype_key = g._genotype_key 
    AND g._strain_key = gs._strain_key 
  ;
