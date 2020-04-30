/*
 View "mgd.all_allele_view"
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
    a.ismixed,   -- gender
    a._refs_key,
    a._markerallele_status_key,
    m.symbol AS markersymbol,
    t.term,
    t.sequencenum AS statusnum,
    s.strain,
    t2.term AS collection,
    t3.term AS markerallele_status,
    r.numericpart AS jnum,
    r.jnumid,
    r.citation,
    r.short_citation
   FROM all_allele a
     LEFT JOIN mrk_marker m ON a._marker_key = m._marker_key
     LEFT JOIN bib_citation_cache r ON a._refs_key = r._refs_key
     join prb_strain s on a._strain_key = s._strain_key
     join voc_term t on a._allele_status_key = t._term_key 
     join voc_term t2 on a._collection_key = t2._term_key 
     join voc_term t3 on a._markerallele_status_key = t3._term_key

  ;

