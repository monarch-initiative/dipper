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
*/
 SELECT e._annotevidence_key,
    e._annot_key,
    e._evidenceterm_key,
    e._refs_key,
    e.inferredfrom,
    e._createdby_key,
    e._modifiedby_key,
    e.creation_date,
    e.modification_date,
    t.abbreviation AS evidencecode,
    t.sequencenum AS evidenceseqnum,
    c.jnumid,
    c.numericpart AS jnum,
    c.short_citation,
    u1.login AS createdby,
    u2.login AS modifiedby
   FROM voc_evidence e,
    voc_term t,
    bib_citation_cache c,
    mgi_user u1,
    mgi_user u2
  WHERE e._evidenceterm_key = t._term_key 
    AND e._refs_key = c._refs_key
    AND e._createdby_key = u1._user_key 
    AND e._modifiedby_key = u2._user_key
 ;
