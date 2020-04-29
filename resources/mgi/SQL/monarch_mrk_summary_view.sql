/*
mgd=> \d+ mrk_summary_view
                                        View "mgd.mrk_summary_view"
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
 SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.private,
    a.preferred,
    a2.accid AS mgiid,
    mt.name AS subtype,
    (((m.symbol || ', '::text) || m.name) || ', Chr '::text) || m.chromosome AS description,
    m.symbol AS short_description
   FROM acc_accession a,
    acc_accession a2,
    mrk_marker m,
    mrk_types mt
  WHERE m._organism_key = 1 
    AND m._marker_type_key = mt._marker_type_key 
    AND m._marker_key = a._object_key 
    AND a._mgitype_key = 2 
    AND a.private = 0 
    AND a._object_key = a2._object_key 
    AND a2._logicaldb_key = 1 
    AND a2._mgitype_key = 2 
    AND a2.prefixpart = 'MGI:'::text 
    AND a2.preferred = 1
  ;
