/*
mgd=> \d+ all_allele_view
                                            View "mgd.all_allele_view"
          Column          |            Type             | Collation | Nullable | Default | Storage  | Description 
--------------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _allele_key              | integer                     |           |          |         | plain    | 
 _marker_key              | integer                     |           |          |         | plain    | 
 _strain_key              | integer                     |           |          |         | plain    | 
 _mode_key                | integer                     |           |          |         | plain    | 
 _allele_type_key         | integer                     |           |          |         | plain    | 
 _allele_status_key       | integer                     |           |          |         | plain    | 
 _transmission_key        | integer                     |           |          |         | plain    | 
 _collection_key          | integer                     |           |          |         | plain    | 
 symbol                   | text                        |           |          |         | extended | 
 name                     | text                        |           |          |         | extended | 
 iswildtype               | smallint                    |           |          |         | plain    | 
 isextinct                | smallint                    |           |          |         | plain    | 
 ismixed                  | smallint                    |           |          |         | plain    | 
 _refs_key                | integer                     |           |          |         | plain    | 
 _markerallele_status_key | integer                     |           |          |         | plain    | 
 _createdby_key           | integer                     |           |          |         | plain    | 
 _modifiedby_key          | integer                     |           |          |         | plain    | 
 _approvedby_key          | integer                     |           |          |         | plain    | 
 approval_date            | timestamp without time zone |           |          |         | plain    | 
 creation_date            | timestamp without time zone |           |          |         | plain    | 
 modification_date        | timestamp without time zone |           |          |         | plain    | 
 markersymbol             | text                        |           |          |         | extended | 
 term                     | text                        |           |          |         | extended | 
 statusnum                | integer                     |           |          |         | plain    | 
 strain                   | text                        |           |          |         | extended | 
 collection               | text                        |           |          |         | extended | 
 createdby                | text                        |           |          |         | extended | 
 modifiedby               | text                        |           |          |         | extended | 
 approvedby               | text                        |           |          |         | extended | 
 markerallele_status      | text                        |           |          |         | extended | 
 jnum                     | integer                     |           |          |         | plain    | 
 jnumid                   | text                        |           |          |         | extended | 
 citation                 | text                        |           |          |         | extended | 
 short_citation           | text                        |           |          |         | extended | 
View definition:
*/

 SELECT a._allele_key,
    a._marker_key,
    a._strain_key,
    a._mode_key,
    a._allele_type_key,
    a._allele_status_key,
    a._transmission_key,
    a._collection_key,
    a.symbol,
    a.name,
    a.iswildtype,
    a.isextinct,
    a.ismixed,
    a._refs_key,
    a._markerallele_status_key,
    a._createdby_key,
    a._modifiedby_key,
    a._approvedby_key,
    a.approval_date,
    a.creation_date,
    a.modification_date,
    m.symbol AS markersymbol,
    t.term,
    t.sequencenum AS statusnum,
    s.strain,
    t2.term AS collection,
    u1.login AS createdby,
    u2.login AS modifiedby,
    u3.login AS approvedby,
    t3.term AS markerallele_status,
    r.numericpart AS jnum,
    r.jnumid,
    r.citation,
    r.short_citation
   FROM all_allele a
     LEFT JOIN mrk_marker m ON a._marker_key = m._marker_key
     LEFT JOIN bib_citation_cache r ON a._refs_key = r._refs_key
     LEFT JOIN mgi_user u3 ON a._approvedby_key = u3._user_key,
    prb_strain s,
    voc_term t,
    voc_term t2,
    voc_term t3,
    mgi_user u1,
    mgi_user u2
  WHERE a._allele_status_key = t._term_key 
    AND a._collection_key = t2._term_key 
    AND a._markerallele_status_key = t3._term_key 
    AND a._strain_key = s._strain_key 
    AND a._createdby_key = u1._user_key 
    AND a._modifiedby_key = u2._user_key
  ;

