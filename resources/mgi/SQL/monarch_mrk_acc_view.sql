/*
First derived from View "mgd.mrk_acc_view"
~ 2 minutes
*/


SELECT
	mav._accession_key,
    mav.accid,
    mav.prefixpart,
    mav.numericpart,
    --a._logicaldb_key,
    mav._object_key as mgi_internal,
    --a._mgitype_key,
    at.name as mgitype,
    --a.private,
    --a.preferred,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,
    mav.logicaldb,
    --m._organism_key
    mov.latinname

FROM mrk_acc_view mav
join acc_mgitype at on mav._mgitype_key = at._mgitype_key
join mgi_organism_view mov on mav._organism_key =  mov._organism_key

WHERE mav.private !=1
  AND mav.preferred = 1
;
