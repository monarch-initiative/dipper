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
 SELECT 
    v._annot_key,
    acc_type.name,
    q.abbreviation AS qualifier,
    t0.term,
    t0._logicaldb_key, -- voc_term_view
    t1.term as annottype_term

    a._mgitype_key,
    a._evidencevocab_key,
    a.name AS annottype
   FROM voc_annot v,
    join acc_type on  v._object_key = _mgitype_key._mgitype_key
    join voc_term_view t0 on v._term_key = t0._term_key 
    join voc_annottype a on v._annottype_key = a._annottype_key 
    join voc_term q on v._qualifier_key = q._term_key
    join voc_term_view t1 on a._vocab_key = t1._term_key
  WHERE 
    AND 
    AND 
  ;
