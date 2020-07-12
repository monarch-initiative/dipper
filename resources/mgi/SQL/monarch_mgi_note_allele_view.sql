/*
First derived from View "mgd.mgi_note_allele_view"
 ~ 2 ninutes
*/

 SELECT
	nav._note_key,
    nav._object_key as mgi_internal,
    --n._mgitype_key,
    --n._notetype_key,
    --n._createdby_key,
    --n._modifiedby_key,
    --n.creation_date,
    --n.modification_date,
    nav.notetype,
    nav.note,
    --c.sequencenum,
    nav.mgitype
    --u1.login AS createdby,
    --u2.login AS modifiedby
FROM mgi_note_allele_view nav
;
