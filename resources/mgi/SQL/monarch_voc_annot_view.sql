/*
First derived from View "mgd.voc_annot_view"

~ 45 seconds   -> killed in 11 minutes
*/

 SELECT
 	vav._annot_key,
    -- vav._annottype_key,
    at.name as mgitype,  -- vav._mgitype_key,
    vav._object_key as mgi_internal,  

    --vav._term_key,
    vav.term, -- voc_term_view._term_key/term,

    --vav._qualifier_key,
    vav.qualifier, -- voc_term._term_key/abbreviation ,

    --vav._vocab_key,
    vav.annottype, --a.name

    --vav._evidencevocab_key,
    vvv.name as evidence

FROM voc_annot_view vav
	join acc_mgitype at on vav._mgitype_key = at._mgitype_key
	join voc_vocab_view vvv on vav._evidencevocab_key = vvv._vocab_key
  ;
