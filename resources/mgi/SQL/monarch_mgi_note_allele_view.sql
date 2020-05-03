/*
First derived from View "mgd.mgi_note_allele_view"
 ~ 35 seconds
*/

 SELECT n._note_key,
    at.name || ':' || n._object_key as mgi_internal,
    t.notetype,
    c.note,
    c.sequencenum,
    m.name AS mgitype

   FROM mgi_note n
    join mgi_notetype t on n._notetype_key = t._notetype_key 
    join mgi_notechunk c on n._note_key = c._note_key 
    join acc_mgitype m on n._mgitype_key = m._mgitype_key 
    join acc_mgitype at on t._mgitype_key =  at._mgitype_key
  WHERE at.name = 'Allele'
    AND t.private != 1

 ;
