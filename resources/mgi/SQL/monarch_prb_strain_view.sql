/*
View "mgd.prb_strain_view"
*/
 SELECT s._strain_key,
    s._species_key,
    s._straintype_key,
    s.strain,
    s.standard,
    s.geneticbackground,
    sp.term AS species,
    st.term AS straintype,
   FROM prb_strain s
    join voc_term sp on s._species_key = sp._term_key 
    join voc_term st on s._straintype_key = st._term_key 
  WHERE s.private = 0
  ;


