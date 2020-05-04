/*
First derived from View "mgd.gxd_genotype_view"

~ 5 seconds
*/

 SELECT
	ggv._genotype_key,
    --ggv._strain_key,
    ggv.isconditional,
    ggv.note,
    --ggv._existsas_key,
    --g._createdby_key,
    --g._modifiedby_key,
    --g.creation_date,
    --g.modification_date,
    ggv.strain,
    ggv.mgiid,
    ggv.displayit,
    --u1.login AS createdby,
    --u2.login AS modifiedby,
    ggv.existsas
   FROM gxd_genotype_view ggv
;
