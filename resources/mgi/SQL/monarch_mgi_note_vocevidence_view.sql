/*
 First derived from View "mgd.mgi_note_vocevidence_view"

 ~ 20 seconds
*/

 SELECT n._note_key,
    'AnnotationEvidence:' || n._object_key as mgi_internal,
    n._mgitype_key,
    n._notetype_key,
    t.notetype,
    c.note,
    c.sequencenum,
    at.name AS mgitype
   FROM mgi_note n
    join mgi_notetype t on n._notetype_key = t._notetype_key 
    join mgi_notechunk c on n._note_key = c._note_key 
    join acc_mgitype at on  n._mgitype_key = at._mgitype_key 

  WHERE at.name = 'Annotation Evidence'
 ;

