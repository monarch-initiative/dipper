/*
mgd=> \d+ voc_evidence_view
                                        View "mgd.voc_evidence_view"
       Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
--------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _annotevidence_key | integer                     |           |          |         | plain    | 
 _annot_key         | integer                     |           |          |         | plain    | 
 _evidenceterm_key  | integer                     |           |          |         | plain    | 
 _refs_key          | integer                     |           |          |         | plain    | 
 inferredfrom       | text                        |           |          |         | extended | 
 _createdby_key     | integer                     |           |          |         | plain    | 
 _modifiedby_key    | integer                     |           |          |         | plain    | 
 creation_date      | timestamp without time zone |           |          |         | plain    | 
 modification_date  | timestamp without time zone |           |          |         | plain    | 
 evidencecode       | text                        |           |          |         | extended | 
 evidenceseqnum     | integer                     |           |          |         | plain    | 
 jnumid             | text                        |           |          |         | extended | 
 jnum               | integer                     |           |          |         | plain    | 
 short_citation     | text                        |           |          |         | extended | 
 createdby          | text                        |           |          |         | extended | 
 modifiedby         | text                        |           |          |         | extended | 
View definition:

TEC: unused?

*/
 SELECT e._annotevidence_key,
    e._annot_key,
    e.inferredfrom,
    t0.abbreviation AS evidencecode,
    t0.sequencenum AS evidenceseqnum,
    t1.term as bibreftyp,
    c.jnumid,
    c.numericpart AS jnum,
    c.mgiid as as mgipub,
    c.pubmedid,
    c.doiid,
    c.short_citation
   FROM voc_evidence e
    join voc_term t0 on e._evidenceterm_key = t0._term_key
    join bibs_ref on  e._refs_key = bibs_ref._refs_key 
    join voc_term t1 on bibs_ref._refferencetype_key = t1._term_key
    join bib_citation_cache c on e._refs_key = c._refs_key
 ;
