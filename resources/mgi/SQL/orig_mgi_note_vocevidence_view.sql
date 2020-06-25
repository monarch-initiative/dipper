/*
mgd=> \d+ mgi_note_vocevidence_view
                                   View "mgd.mgi_note_vocevidence_view"
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
    n._createdby_key,
    n._modifiedby_key,
    n.creation_date,
    n.modification_date,
    t.notetype,
    c.note,
    c.sequencenum,
    m.name AS mgitype,
    u1.login AS createdby,
    u2.login AS modifiedby
   FROM mgi_note n,
    mgi_notetype t,
    mgi_notechunk c,
    acc_mgitype m,
    mgi_user u1,
    mgi_user u2
  WHERE n._notetype_key = t._notetype_key 
    AND t._mgitype_key = 25 
    AND n._note_key = c._note_key 
    AND n._mgitype_key = m._mgitype_key 
    AND n._createdby_key = u1._user_key 
    AND n._modifiedby_key = u2._user_key
 ;

