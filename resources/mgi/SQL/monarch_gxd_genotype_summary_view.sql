/*
View "mgd.gxd_genotype_summary_view"
*/
 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.private,
    a.preferred,
    a.accid AS mgiid,
    s.strain AS subtype,
    (((s.strain || ' '::text) || a1.symbol) || ','::text) || a2.symbol AS description,
    (a1.symbol || ','::text) || a2.symbol AS short_description,
    l.name AS logicaldb
   FROM acc_accession a
     JOIN gxd_genotype g ON a._object_key = g._genotype_key
     JOIN prb_strain s ON g._strain_key = s._strain_key
     JOIN gxd_allelepair ap ON g._genotype_key = ap._genotype_key
     JOIN all_allele a1 ON ap._allele_key_1 = a1._allele_key
     JOIN all_allele a2 ON ap._allele_key_2 = a2._allele_key
     JOIN acc_logicaldb l ON a._logicaldb_key = l._logicaldb_key
  WHERE a._mgitype_key = 12 AND a._logicaldb_key = 1 AND a.prefixpart = 'MGI:'::text AND a.preferred = 1 AND a._logicaldb_key = l._logicaldb_key
UNION
 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.private,
    a.preferred,
    a.accid AS mgiid,
    s.strain AS subtype,
    (s.strain || ' '::text) || a1.symbol AS description,
    a1.symbol AS short_description,
    l.name AS logicaldb
   FROM acc_accession a
     JOIN gxd_genotype g ON a._object_key = g._genotype_key
     JOIN prb_strain s ON g._strain_key = s._strain_key
     JOIN gxd_allelepair ap ON g._genotype_key = ap._genotype_key AND ap._allele_key_2 IS NULL
     JOIN all_allele a1 ON ap._allele_key_1 = a1._allele_key
     JOIN acc_logicaldb l ON a._logicaldb_key = l._logicaldb_key
  WHERE a._mgitype_key = 12 AND a._logicaldb_key = 1 AND a.prefixpart = 'MGI:'::text AND a.preferred = 1 AND a._logicaldb_key = l._logicaldb_key
UNION
 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.private,
    a.preferred,
    a.accid AS mgiid,
    s.strain AS subtype,
    s.strain AS description,
    NULL::text AS short_description,
    l.name AS logicaldb
   FROM acc_accession a
     JOIN gxd_genotype g ON a._object_key = g._genotype_key
     JOIN prb_strain s ON g._strain_key = s._strain_key
     JOIN acc_logicaldb l ON a._logicaldb_key = l._logicaldb_key
  WHERE a._mgitype_key = 12 
    AND a._logicaldb_key = 1 
    AND a.prefixpart = 'MGI:'::text 
    AND a.preferred = 1 
    AND a._logicaldb_key = l._logicaldb_key 
    AND NOT (EXISTS ( 
        SELECT 1
        FROM gxd_allelepair ap
        WHERE g._genotype_key = ap._genotype_key))
;
