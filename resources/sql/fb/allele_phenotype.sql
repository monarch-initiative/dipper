/* Query to fetch allele to phenotype associations

   This is similar to the file here:
   ftp://ftp.flybase.net/releases/current/precomputed_files/
   alleles/allele_phenotypic_data_fb_2019_02.tsv.gz

   Instead we're opting to parse the text in derived_pheno_class
   because it contains IDs for FBcv terms and alleles

   Note that the above file is larger because it also includes
   cvterm.name = 'derived_pheno_manifest'
   which contain anatomy terms
 */
SELECT
    feature.uniquename as allele_id,
    featureprop.value as pheno_desc,
    cvterm.name as pheno_type,
    pub.uniquename as pub_id,
    pub.title as pub_title,
    pmid.accession as pmid_id

FROM feature

JOIN featureprop
ON feature.feature_id = featureprop.feature_id

JOIN cvterm
ON featureprop.type_id = cvterm.cvterm_id
    AND (cvterm.name = 'derived_pheno_class' OR cvterm.name = 'derived_pheno_manifest')

JOIN featureprop_pub
ON featureprop.featureprop_id = featureprop_pub.featureprop_id

JOIN pub
ON featureprop_pub.pub_id = pub.pub_id

LEFT OUTER JOIN (

    SELECT pub_lookup.pub_id, xref.db_id, xref.dbxref_id, xref.accession, dbinfo.name
    FROM pub_dbxref as pub_lookup

    JOIN dbxref xref
    ON pub_lookup.dbxref_id = xref.dbxref_id

    JOIN db dbinfo
    ON xref.db_id = dbinfo.db_id AND dbinfo.name = 'pubmed'

) pmid ON pub.pub_id = pmid.pub_id
