/*
View "mgd.mrk_marker_view"
*/
 SELECT m._marker_key,
    m._organism_key,
    m._marker_status_key,
    m._marker_type_key,
    m.symbol,
    m.name,
    m.chromosome,
    m.cytogeneticoffset,
    m.cmoffset,
    ((s.commonname || ' ('::text) || s.latinname) || ')'::text AS organism,
    s.commonname,
    s.latinname,
    ms.status,
    mt.name AS markertype

   FROM mrk_marker m 
    join mgi_organism s on m._organism_key = s._organism_key 
    join mrk_status ms on m._marker_status_key = ms._marker_status_key 
    join mrk_types mt on m._marker_type_key = mt._marker_type_key

  ;
