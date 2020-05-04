/*
First derived from View "mgd.prb_strain_view"

~ 5 seconds  -> another sloooow one ...

TODO retest in case they just have an overloaded werver

*/

SELECT
	psv._strain_key,
    --s._species_key,
    --s._straintype_key,
    psv.strain,
    psv.standard,
    --s.private,
    psv.geneticbackground,
    --s._createdby_key,
    --s._modifiedby_key,
    --s.creation_date,
    --s.modification_date,
    psv.species,
    psv.straintype
    --u1.login AS createdby,
    --u2.login AS modifiedby
FROM prb_strain_view psv
WHERE psv.private != 1
;

/*
 SELECT 
 	s._strain_key,
    s.strain,
    s.standard,
    s.geneticbackground,
    sp.term AS strain_species,
    st.term AS strain_type
   FROM prb_strain s
    join voc_term sp on s._species_key = sp._term_key 
    join voc_vocab spv on sp._vocab_key = spv._vocab_key 
    join voc_term st on s._straintype_key = st._term_key 
    join voc_vocab stv on sp._vocab_key = stv._vocab_key 
  WHERE s.private != 1
  ;
*/