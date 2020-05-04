/*
	~45 seconds
*/

select 
 msv._accession_key, 
 msv.accid,
 msv.prefixpart, 
 msv.numericpart,
 al.name as logicaldb,
 at.name as mgi_type,
 msv._object_key as mgi_internal,
 msv.preferred,
 --_createdby_key 
 --_modifiedby_key
 --creation_date
 --modification_date
 msv.mgiid,  
 msv.subtype, 
 msv.description, 
 msv.short_description

from  mrk_summary_view msv
 join acc_mgitype at on msv._mgitype_key = at._mgitype_key
 join acc_logicaldb_view al on msv._logicaldb_key = al._logicaldb_key

where msv.private != 1
  and msv.preferred = 1
;