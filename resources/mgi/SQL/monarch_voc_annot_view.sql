/*
View "mgd.voc_annot_view"
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
  ;
