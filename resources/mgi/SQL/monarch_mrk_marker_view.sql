/*
mgd=> \d+ mrk_marker_view
                                         View "mgd.mrk_marker_view"
       Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
--------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _marker_key        | integer                     |           |          |         | plain    | 
 _organism_key      | integer                     |           |          |         | plain    | 
 _marker_status_key | integer                     |           |          |         | plain    | 
 _marker_type_key   | integer                     |           |          |         | plain    | 
 symbol             | text                        |           |          |         | extended | 
 name               | text                        |           |          |         | extended | 
 chromosome         | text                        |           |          |         | extended | 
 cytogeneticoffset  | text                        |           |          |         | extended | 
 cmoffset           | double precision            |           |          |         | plain    | 
 _createdby_key     | integer                     |           |          |         | plain    | 
 _modifiedby_key    | integer                     |           |          |         | plain    | 
 creation_date      | timestamp without time zone |           |          |         | plain    | 
 modification_date  | timestamp without time zone |           |          |         | plain    | 
 organism           | text                        |           |          |         | extended | 
 commonname         | text                        |           |          |         | extended | 
 latinname          | text                        |           |          |         | extended | 
 status             | text                        |           |          |         | extended | 
 markertype         | text                        |           |          |         | extended | 
 createdby          | text                        |           |          |         | extended | 
 modifiedby         | text                        |           |          |         | extended | 
View definition:
*/
 SELECT m._marker_key,
    m._organism_key,
    m._marker_status_key,
    m._marker_type_key,
    m.symbol,
    m.name,
    m.chromosome,
    m.cytogeneticoffset,
    m.cmoffset,
    ((s.commonname || ' ('::text) || s.latinname) || ')'::text AS organism,
    s.commonname,
    s.latinname,
    ms.status,
    mt.name AS markertype

   FROM mrk_marker m,
    mgi_organism s,
    mrk_status ms,
    mrk_types mt

  WHERE m._organism_key = s._organism_key 
  	AND m._marker_status_key = ms._marker_status_key 
  	AND m._marker_type_key = mt._marker_type_key
  ;
