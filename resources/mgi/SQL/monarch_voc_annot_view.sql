/*
View "mgd.voc_annot_view"
*/
 SELECT 
    v._annot_key,
    at.name as annot_type,
    q.abbreviation AS qualifier,
    t0.term annot_name,
    t0._logicaldb_key, -- voc_term_view
    t1.term as annottype_term,
    ev.name as evidence,
    a.name AS annottype
   FROM voc_annot v,
    join acc_type at on v._object_key = at._mgitype_key
    join voc_term_view t0 on v._term_key = t0._term_key 
    join voc_annottype a on v._annottype_key = a._annottype_key 
    join voc_term q on v._qualifier_key = q._term_key
    join voc_term_view t1 on a._vocab_key = t1._term_key 
    join acc_mgitype at on a._mgitype_key =  at._mgitype_key
    join voc_vocab ev on a._evidencevocab_key = ev._vocab_key
  WHERE 
  ;
