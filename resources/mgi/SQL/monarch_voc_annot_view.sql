/*
First derived from View "mgd.voc_annot_view"

~ 45 seconds
*/

 SELECT 
 	vav._annot_key,
    -- vav._annottype_key,
    at.name as mgi_type,  -- vav._mgitype_key,
    vav._object_key as mgi_internal,  

    --vav._term_key,
    vav.term, -- voc_term_view._term_key/term,

    --vav._qualifier_key,
    vav.qualifier, -- voc_term._term_key/abbreviation ,

    --vav._vocab_key,
    vav.annottype, --a.name

    --vav._evidencevocab_key,
    aev.name as evidence_type

FROM voc_annot_view vav
	join acc_mgitype at on vav._mgitype_key = at._mgitype_key
	join voc_vocab aev on vav._evidencevocab_key = aev._vocab_key

  ;
