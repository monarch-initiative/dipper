/*
mgd=> \d+ prb_strain_genotype_view
                                     View "mgd.prb_strain_genotype_view"
       Column        |            Type             | Collation | Nullable | Default | Storage  | Description 
---------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _straingenotype_key | integer                     |           |          |         | plain    | 
 _strain_key         | integer                     |           |          |         | plain    | 
 _genotype_key       | integer                     |           |          |         | plain    | 
 _qualifier_key      | integer                     |           |          |         | plain    | 
 _createdby_key      | integer                     |           |          |         | plain    | 
 _modifiedby_key     | integer                     |           |          |         | plain    | 
 creation_date       | timestamp without time zone |           |          |         | plain    | 
 modification_date   | timestamp without time zone |           |          |         | plain    | 
 strain              | text                        |           |          |         | extended | 
 description         | text                        |           |          |         | extended | 
 mgiid               | text                        |           |          |         | extended | 
 qualifier           | text                        |           |          |         | extended | 
 modifiedby          | text                        |           |          |         | extended | 
View definition:
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
    u.login AS modifiedby
   FROM prb_strain_genotype s,
    prb_strain ss,
    acc_accession a,
    gxd_genotype g,
    prb_strain gs,
    voc_term t,
    mgi_user u
  WHERE s._strain_key = ss._strain_key 
    AND s._qualifier_key = t._term_key 
    AND s._genotype_key = a._object_key 
    AND a._mgitype_key = 12 
    AND a._logicaldb_key = 1 
    AND s._genotype_key = g._genotype_key 
    AND g._strain_key = gs._strain_key 
    AND s._modifiedby_key = u._user_key
  ;
