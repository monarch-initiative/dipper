/*
View "mgd.prb_strain_view"
*/
 SELECT s._strain_key,
    s.strain,
    s.standard,
    s.geneticbackground,
    sp.term AS strain_species,
    st.term AS strain_type,
   FROM prb_strain s
    join voc_term sp on s._species_key = sp._term_key 
    join voc_vocab spv on sp._vocab_key = spv._vocab_key 
    join voc_term st on s._straintype_key = st._term_key 
    join voc_vocab stv on sp._vocab_key = stv._vocab_key 
  WHERE s.private != 1
    AND spv.term = 'Strain Species'
    AND stv.term = 'Strain Type'
  ;
