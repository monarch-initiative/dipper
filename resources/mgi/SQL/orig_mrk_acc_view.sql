/*
mgd=> \d+ mrk_acc_view
                                          View "mgd.mrk_acc_view"
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
 logicaldb         | text                        |           |          |         | extended | 
 _organism_key     | integer                     |           |          |         | plain    | 
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
    a._createdby_key,
    a._modifiedby_key,
    a.creation_date,
    a.modification_date,
    l.name AS logicaldb,
    m._organism_key
   FROM acc_accession a,
    acc_logicaldb l,
    mrk_marker m
  WHERE a._mgitype_key = 2 
    AND a._logicaldb_key = l._logicaldb_key 
    AND a._object_key = m._marker_key
 ;