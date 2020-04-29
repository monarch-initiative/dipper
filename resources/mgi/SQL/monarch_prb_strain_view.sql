/*
mgd=> \d+ prb_strain_view
                                        View "mgd.prb_strain_view"
      Column       |            Type             | Collation | Nullable | Default | Storage  | Description 
-------------------+-----------------------------+-----------+----------+---------+----------+-------------
 _strain_key       | integer                     |           |          |         | plain    | 
 _species_key      | integer                     |           |          |         | plain    | 
 _straintype_key   | integer                     |           |          |         | plain    | 
 strain            | text                        |           |          |         | extended | 
 standard          | smallint                    |           |          |         | plain    | 
 private           | smallint                    |           |          |         | plain    | 
 geneticbackground | smallint                    |           |          |         | plain    | 
 _createdby_key    | integer                     |           |          |         | plain    | 
 _modifiedby_key   | integer                     |           |          |         | plain    | 
 creation_date     | timestamp without time zone |           |          |         | plain    | 
 modification_date | timestamp without time zone |           |          |         | plain    | 
 species           | text                        |           |          |         | extended | 
 straintype        | text                        |           |          |         | extended | 
 createdby         | text                        |           |          |         | extended | 
 modifiedby        | text                        |           |          |         | extended | 
View definition:
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


