/*
First derived from View "mgd.gxd_allelepair_view"

 ~ 1 minute

*/

SELECT
	gav._allelepair_key,
    --gav._genotype_key,
    ggv.strain,
    ggv.mgiid,
    --gav._allele_key_1,
    --gav._allele_key_2,
    --gav._marker_key,
    --gav._mutantcellline_key_1,
    cl1.cellline as mutantcellline1,
    --gav._mutantcellline_key_2,
    cl2.cellline as mutantcellline2,
    --gav._pairstate_key,
    --gav._compound_key,
    gav.sequencenum,
    --a._createdby_key,
    --a._modifiedby_key,
    --a.creation_date,
    --a.modification_date,

    gav.symbol,
    gav.chromosome,
    gav.allele1,
    gav.allele2,
    gav.allelestate,
    gav.compound

FROM gxd_allelepair_view gav
join gxd_genotype_view ggv on gav._genotype_key = ggv._genotype_key
left join all_cellline_view cl1 on gav._mutantcellline_key_1 = cl1._cellline_key
left join all_cellline_view cl2 on gav._mutantcellline_key_2 = cl2._cellline_key
;
