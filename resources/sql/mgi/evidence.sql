SELECT
  voc_evidence_view._annotevidence_key,
  voc_evidence_view._annot_key,
  voc_evidence_view.evidencecode,
  voc_evidence_view.jnumid,
  voc_term.term,
  voc_evidence_property.value,
  voc_annot_view.annottype
FROM voc_evidence_view
JOIN voc_annot_view
  ON voc_evidence_view._annot_key = voc_annot_view._annot_key
LEFT OUTER JOIN voc_evidence_property
  ON voc_evidence_view._annotevidence_key = voc_evidence_property._annotevidence_key
LEFT OUTER JOIN voc_term
  ON voc_evidence_property._propertyterm_key = voc_term._term_key
WHERE voc_annot_view.annottype <> 'GO/Marker'