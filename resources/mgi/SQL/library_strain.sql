-- http://www.informatics.jax.org/schema_pg/tables/prb_probe.html
--      _segmenttype_key?
-- http://www.informatics.jax.org/schema_pg/tables/prb_source.html
-- http://www.informatics.jax.org/schema_pg/tables/prb_strain.html

select distinct -- count(*)
    --
    latinname,
    --
    strain,
    --
    strain_type.term strain_type,
    --
    gender.term gender,
    --
    tissue,
    --
    geneticbackground,
    --
    prb_refs.term prb_refs
    
from prb_source  -- library?
    -- strain label
    join prb_strain on prb_source._strain_key = prb_strain._strain_key

    -- species indicator for taxon
    join mgi_organism on prb_source._organism_key = mgi_organism._organism_key

    -- strain type-hint 
    join voc_term strain_type on prb_source._straintype_key = strain_type._term_key

    -- gender
    join voc_term gender on prb_source._gender_key = gender._term_key
    
    -- anatomy-ish
    left join prb_tissue on prb_source._tissue_key = prb_tissue._tissue_key 

    -- probe_library_publications?
    join bib_refs on prb_source._refs_key = bib_refs._refs_key
    -- via voc_term
    left join voc_term prb_refs on bib_refs._refs_key = prb_refs._term_key

 where prb_strain.standard = 1
   and prb_strain.private != 1 
   and prb_tissue.standard = 1  
--
limit 100
;
