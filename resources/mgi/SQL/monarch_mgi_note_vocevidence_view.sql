/*
 First derived from View "mgd.mgi_note_vocevidence_view"

 ~ 30 seconds
*/

 SELECT
	nvv._note_key,
    nvv._object_key as mgi_internal,
    --n._mgitype_key,
    --n._notetype_key,
    --n._createdby_key,
    --n._modifiedby_key,
    --n.creation_date,
    --n.modification_date,
    nvv.notetype,
    nvv.note,
    --c.sequencenum,
    nvv.mgitype
    --u1.login AS createdby,
    --u2.login AS modifiedby

   FROM mgi_note_vocevidence_view nvv
;
