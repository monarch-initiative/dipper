SELECT feature_relationship_id, subject_id, object_id, name, rank, value
FROM feature_relationship
INNER JOIN cvterm
ON feature_relationship.type_id = cvterm.cvterm_id