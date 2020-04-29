/*
mgd=> \d+ mgi_note_allele_view
                                      View "mgd.mgi_note_allele_view"
      Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
-------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _note_key         | integer                     |           |          |         | plain    | 
 _object_key       | integer                     |           |          |         | plain    | 
 _mgitype_key      | integer                     |           |          |         | plain    | 
 _notetype_key     | integer                     |           |          |         | plain    | 
 _createdby_key    | integer                     |           |          |         | plain    | 
 _modifiedby_key   | integer                     |           |          |         | plain    | 
 creation_date     | timestamp without time zone |           |          |         | plain    | 
 modification_date | timestamp without time zone |           |          |         | plain    | 
 notetype          | text                        |           |          |         | extended | 
 note              | text                        |           |          |         | extended | 
 sequencenum       | integer                     |           |          |         | plain    | 
 mgitype           | text                        |           |          |         | extended | 
 createdby         | text                        |           |          |         | extended | 
 modifiedby        | text                        |           |          |         | extended | 
View definition:
*/

 SELECT n._note_key,
    n._object_key,
    n._mgitype_key,
    n._notetype_key,
    t.notetype,
    c.note,
    c.sequencenum,
    m.name AS mgitype

   FROM mgi_note n,
    mgi_notetype t,
    mgi_notechunk c,
    acc_mgitype m

  WHERE n._notetype_key = t._notetype_key 
    AND t._mgitype_key = 11 
    AND n._note_key = c._note_key 
    AND n._mgitype_key = m._mgitype_key 
 ;
