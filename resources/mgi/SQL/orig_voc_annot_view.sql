/*
mgd=> \d+ voc_annot_view
                                         View "mgd.voc_annot_view"
       Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
--------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _annot_key         | integer                     |           |          |         | plain    | 
 _annottype_key     | integer                     |           |          |         | plain    | 
 _object_key        | integer                     |           |          |         | plain    | 
 _term_key          | integer                     |           |          |         | plain    | 
 _qualifier_key     | integer                     |           |          |         | plain    | 
 creation_date      | timestamp without time zone |           |          |         | plain    | 
 modification_date  | timestamp without time zone |           |          |         | plain    | 
 qualifier          | text                        |           |          |         | extended | 
 term               | text                        |           |          |         | extended | 
 sequencenum        | integer                     |           |          |         | plain    | 
 accid              | text                        |           |          |         | extended | 
 _logicaldb_key     | integer                     |           |          |         | plain    | 
 _vocab_key         | integer                     |           |          |         | plain    | 
 _mgitype_key       | integer                     |           |          |         | plain    | 
 _evidencevocab_key | integer                     |           |          |         | plain    | 
 annottype          | text                        |           |          |         | extended | 
View definition:
*/
 SELECT v._annot_key,
    v._annottype_key,
    v._object_key,
    v._term_key,
    v._qualifier_key,
    v.creation_date,
    v.modification_date,
    q.abbreviation AS qualifier,
    t.term,
    t.sequencenum,
    t.accid,
    t._logicaldb_key,
    a._vocab_key,
    a._mgitype_key,
    a._evidencevocab_key,
    a.name AS annottype
   FROM voc_annot v,
    voc_term_view t,
    voc_annottype a,
    voc_term q
  WHERE v._term_key = t._term_key 
    AND v._annottype_key = a._annottype_key 
    AND v._qualifier_key = q._term_key
  ;
