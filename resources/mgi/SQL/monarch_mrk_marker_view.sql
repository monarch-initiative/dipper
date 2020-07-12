/*
First derived from View "mgd.mrk_marker_view"

~2.5 minutes

*/

SELECT
	mmv._marker_key,
    --m._organism_key,
    --m._marker_status_key,
    --m._marker_type_key,
    mmv.symbol,
    mmv.name,
    mmv.chromosome,
    mmv.cytogeneticoffset,
    mmv.cmoffset,
    --m._createdby_key,
    --m._modifiedby_key,
    --m.creation_date,
    --m.modification_date,
    mmv.organism,
    --mmv.commonname,
    mmv.latinname,
    mmv.status,
    mmv.markertype
    --u1.login AS createdby,
    --u2.login AS modifiedby
FROM mrk_marker_view mmv
;
