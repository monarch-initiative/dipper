/*
                                      View "mgd.all_summary_view"
      Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
-------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _accession_key    | integer                     |           |          |         | plain    | 
 accid             | text                        |           |          |         | extended | 
 prefixpart        | text                        |           |          |         | extended | 
 numericpart       | integer                     |           |          |         | plain    | 
 _logicaldb_key    | integer                     |           |          |         | plain    | 
 _object_key       | integer                     |           |          |         | plain    | 
 _mgitype_key      | integer                     |           |          |         | plain    | 
 private           | smallint                    |           |          |         | plain    | 
 preferred         | smallint                    |           |          |         | plain    | 
 _createdby_key    | integer                     |           |          |         | plain    | 
 _modifiedby_key   | integer                     |           |          |         | plain    | 
 creation_date     | timestamp without time zone |           |          |         | plain    | 
 modification_date | timestamp without time zone |           |          |         | plain    | 
 mgiid             | text                        |           |          |         | extended | 
 subtype           | text                        |           |          |         | extended | 
 description       | text                        |           |          |         | extended | 
 short_description | text                        |           |          |         | extended | 
View definition:
*/

 SELECT 
    'Allele:'::text || al._allele_key as _mgiobjid, -- grouping
    a.accid,
    a.prefixpart,
    a.numericpart,
    a2.accid AS mgiid,

    m._marker_key, 
    m.symbol as refmrkrsym, 
    m.name as refmrkenm, -- refference, or NULL if variant
    s.strain,

    aldb.name as db_name,
    org.commonname as orgcommon,
    org.latinname as  orglatin, 
    al.symbol,  -- AS short_description
    al.name,    -- "symbol, name"  AS description,

    t0.term AS subtype,
    t.term as status,
    --t2.term AS collection,
    t3.term AS markerallele_status,
    t4.term AS transmission,

    r.numericpart AS jnum,
    --r.jnumid,
    --r.citation,
    r.short_citation

   FROM all_allele al
   join acc_accession a on  a._object_key = al._allele_key 
   join acc_accession a2 on a._object_key = a2._object_key 
   join acc_mgitype typ on a._mgitype_key = typ._mgitype_key
   join acc_logicaldb aldb on a._logicaldb_key = aldb._logicaldb_key
   join mgi_organism org on aldb._organism_key  = org._organism_key
   join voc_term t0 on al._allele_type_key = t0._term_key

   -----------------------------------------------------------

   LEFT JOIN mrk_marker m ON al._marker_key = m._marker_key
   LEFT JOIN bib_citation_cache r ON al._refs_key = r._refs_key 
   join prb_strain s on al._strain_key = s._strain_key

   join voc_term t on al._allele_status_key = t._term_key 
   --join voc_term t2 on al._collection_key = t2._term_key 
   join voc_term t3 on al._markerallele_status_key = t3._term_key 
   join voc_term t4 on al._transmission_key = t4._term_key

  WHERE typ.name = 'Allele'  -- a._mgitype_key = 11 
    AND a.preferred = 1 
    AND a.private = 0 
    AND a2._logicaldb_key = 1 
    AND a2._mgitype_key = 11 
    AND a2.prefixpart = 'MGI:'::text 
    AND a2.preferred = 1 
    and t.term != 'Deleted'::text

    -----------------------------------------------------

    and  m.symbol != al.symbol
    --
    order by 1 
    --limit 40
;

