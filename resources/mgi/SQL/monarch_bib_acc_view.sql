/*
First derived from View "mgd.bib_acc_view"

~ 1 minute

*/
SELECT
	bav._accession_key,
    bav.accid,
    bav.prefixpart,
    bav.numericpart,
    --bav._logicaldb_key,
    bav._object_key as mgi_internal,
    --bav._mgitype_key,
    bat.name as mgi_type,
    --a.preferred,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,
    bav.logicaldb

FROM bib_acc_view bav
join acc_mgitype bat on bav._mgitype_key =  bat._mgitype_key

WHERE bav.private != 1
  AND bav.preferred = 1
;