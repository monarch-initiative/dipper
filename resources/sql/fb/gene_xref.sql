SELECT feature.uniquename, dbxref.accession, db.name
FROM feature

JOIN cvterm cvt
ON feature.type_id = cvt.cvterm_id
    AND cvt.name = 'gene'

JOIN featureprop
ON feature.feature_id = featureprop.feature_id
    AND feature.is_obsolete = FALSE

JOIN cvterm
ON featureprop.type_id = cvterm.cvterm_id
    AND cvterm.name = 'derived_gene_model_status'
    AND featureprop.value != 'Withdrawn'

JOIN feature_dbxref
ON feature.feature_id = feature_dbxref.feature_id
    AND feature_dbxref.is_current = TRUE

JOIN dbxref
ON feature_dbxref.dbxref_id = dbxref.dbxref_id

JOIN db
ON dbxref.db_id = db.db_id

WHERE db.name = 'EntrezGene' OR db.name = 'HGNC'
