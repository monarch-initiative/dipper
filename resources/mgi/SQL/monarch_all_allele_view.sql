/*
First derived from  View "mgd.all_allele_view"

 ~ 45 seconds
*/

 SELECT 
 	aav._allele_key,
 
    --aav._marker_key, --nullable
    aav.markersymbol,
 
    --aav._strain_key,
    aav.strain,

    --aav._mode_key,
    tin.term as inheritance,

    --aav._allele_type_key,
    tal.term as allele_type,

    --aav._allele_status_key,
    tas.term as allele_status,

    --aav._transmission_key,
    tat.term as transmission,

    --aav._collection_key,
    aav.symbol,
    aav.name,
    aav.iswildtype,
    --aav.isextinct,
    aav.ismixed as mix_gender,
    --aav._refs_key,

    --aav._markerallele_status_key, --nullable
    aav.term AS markerallele_status,

    --aav._createdby_key,
    --aav._modifiedby_key,
    --aav._approvedby_key,
    --aav.approval_date,
    --aav.creation_date,
    --aav.modification_date,
    --m.symbol AS markersymbol,

    aav.term as allele_status,
    --t.sequencenum AS statusnum, -- for internal ordering 

    --aav.collection,
    --u1.login AS createdby,
    --u2.login AS modifiedby,
    --u3.login AS approvedby,
 
    aav.jnum,
    aav.jnumid,
    aav.citation,
    aav.short_citation

FROM all_allele_view aav
  join voc_term tin on aav._mode_key = tin._term_key 
  join voc_term tal on aav._allele_type_key = tal._term_key 
  join voc_term tas on aav._allele_status_key = tas._term_key 
  join voc_term tat on aav._transmission_key = tat._term_key 
where tas.term != 'deleted'

;