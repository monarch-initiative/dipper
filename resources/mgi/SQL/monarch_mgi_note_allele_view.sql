/*
View "mgd.mgi_note_allele_view"
*/

 SELECT n._note_key,
    n._object_key,
    n._mgitype_key,
    n._notetype_key,
    t.notetype,
    c.note,
    c.sequencenum,
    m.name AS mgitype

   FROM mgi_note n
    join mgi_notetype t on n._notetype_key = t._notetype_key 
    join mgi_notechunk c on n._note_key = c._note_key 
    join acc_mgitype m on n._mgitype_key = m._mgitype_key 

  WHERE t._mgitype_key = 11 

 ;
