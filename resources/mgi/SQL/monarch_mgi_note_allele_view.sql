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

   FROM mgi_note n,
    mgi_notetype t,
    mgi_notechunk c,
    acc_mgitype m

  WHERE n._notetype_key = t._notetype_key 
    AND t._mgitype_key = 11 
    AND n._note_key = c._note_key 
    AND n._mgitype_key = m._mgitype_key 
 ;
