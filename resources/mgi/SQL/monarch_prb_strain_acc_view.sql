/*
First derived from View "mgd.prb_strain_acc_view":
~  5 seconds
*/


SELECT
	psav._accession_key,
    psav.accid,
    psav.prefixpart,
    psav.numericpart,
    --a._logicaldb_key,
    psav._object_key as mgi_internal,
    --a._mgitype_key,
    at.name as mgitype,
    --a.private,
    --a.preferred,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,
    psav.logicaldb
FROM prb_strain_acc_view  psav
join acc_mgitype at on psav._mgitype_key = at._mgitype_key

WHERE psav.private != 1
  AND psav.preferred = 1
;
