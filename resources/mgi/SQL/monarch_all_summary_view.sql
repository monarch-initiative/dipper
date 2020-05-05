/*
First derived from View "mgd.all_summary_view"

  ~2 minutes
*/

 SELECT
	asv._accession_key,
    asv.accid,
    asv.prefixpart,
    asv.numericpart,

    --asv._logicaldb_key,
    al.name AS logicaldb,

    asv._object_key as mgi_internal,
    --asv._mgitype_key,
    at.name as mgi_type,

    --a.private,
    --a.preferred,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,

    asv.mgiid,
    asv.subtype,
    asv.description,
    asv.short_description
FROM  all_summary_view asv
  join acc_logicaldb_view al on asv._logicaldb_key = al._logicaldb_key
  join acc_mgitype at on asv._mgitype_key = at._mgitype_key

WHERE asv.private != 1
  AND asv.preferred = 1
;
