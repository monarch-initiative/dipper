/*
View "mgd.prb_strain_view"
*/
 SELECT s._strain_key,
    s._species_key,
    s._straintype_key,
    s.strain,
    s.standard,
    s.private,
    s.geneticbackground,
    sp.term AS species,
    st.term AS straintype,
   FROM prb_strain s,
    voc_term sp,
    voc_term st,
  WHERE s._species_key = sp._term_key 
    AND s._straintype_key = st._term_key 
  ;


