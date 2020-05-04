/*
First derived from  View "mgd.gxd_genotype_summary_view"

	~30 seconds
*/

 SELECT
	ggsv._accession_key,
    ggsv.accid,
    ggsv.prefixpart,
    ggsv.numericpart,
    --a._logicaldb_key,
    ggsv._object_key as mgi_internal,
    --a._mgitype_key,
    at.name as mgi_type,

    --a.private,
    --a.preferred,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,
    ggsv.mgiid,
    ggsv.subtype,   -- strain
    ggsv.description,
    ggsv.short_description,
    ggsv.logicaldb

FROM gxd_genotype_summary_view ggsv
join acc_mgitype at on ggsv._mgitype_key = at._mgitype_key
WHERE ggsv.private !=1
  AND ggsv.preferred = 1
;
