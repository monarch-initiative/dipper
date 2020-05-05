/*
First derived from View "mgd.voc_evidence_view" 


TEC: unused?

*/
 SELECT e._annotevidence_key,
    vat.name as annot_type,
    e.inferredfrom,
    t0.abbreviation AS evidencecode,
    t0.sequencenum AS evidenceseqnum,
    t1.term as bibreftyp,
    c.jnumid,
    c.numericpart AS jnum,
    c.mgiid as mgipub,
    c.pubmedid,
    c.doiid,
    c.short_citation
   FROM voc_evidence e
    join voc_annot va on e._annot_key = va._annot_key 
    join voc_annottype vat on va._annottyp_key = vat._annottype_key
    join voc_term t0 on e._evidenceterm_key = t0._term_key
    join bibs_ref on  e._refs_key = bibs_ref._refs_key 
    join voc_term t1 on bibs_ref._refferencetype_key = t1._term_key
    join bib_citation_cache c on e._refs_key = c._refs_key
 ;
