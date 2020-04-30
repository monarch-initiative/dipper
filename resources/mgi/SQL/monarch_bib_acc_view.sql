/*
View "mgd.bib_acc_view"
*/
SELECT a._accession_key,
    a.accid,
    a.prefixpart,
    a.numericpart,
    a._logicaldb_key,
    a._object_key,
    a._mgitype_key,
    a.preferred,
    l.name AS logicaldb
  FROM acc_accession a 
    join acc_logicaldb l on a._logicaldb_key = l._logicaldb_key
    join acc_mgitype at on a._mgitype_key =  at._mgitype_key
    WHERE at.name = 'Reference'
      AND a.private = 0

;

