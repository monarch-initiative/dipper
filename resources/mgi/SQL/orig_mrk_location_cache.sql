/*
mgd=> \d+ mrk_location_cache
                                                              Table "mgd.mrk_location_cache"
      Column       |            Type             | Collation | Nullable | Default | Storage  | Stats target |                 Description                 
-------------------+-----------------------------+-----------+----------+---------+----------+--------------+---------------------------------------------
 _marker_key       | integer                     |           | not null |         | plain    |              | primary key/foreign key to MRK_Marker
 _marker_type_key  | integer                     |           | not null |         | plain    |              | foreign key to MRK_Types
 _organism_key     | integer                     |           | not null |         | plain    |              | foreign key to MGI_Organism
 chromosome        | text                        |           | not null |         | extended |              | chromosome
 sequencenum       | integer                     |           | not null |         | plain    |              | chromosome order
 cytogeneticoffset | text                        |           |          |         | extended |              | cytogenetic cmOffset
 cmoffset          | double precision            |           |          |         | plain    |              | cytogenetic cmOffset
 genomicchromosome | text                        |           |          |         | extended |              | chromosome
 startcoordinate   | double precision            |           |          |         | plain    |              | start genome coordinate
 endcoordinate     | double precision            |           |          |         | plain    |              | end genome coordinate
 strand            | character(1)                |           |          |         | extended |              | strand for genome coordinate
 mapunits          | text                        |           |          |         | extended |              | units of genome coordinates
 provider          | text                        |           |          |         | extended |              | data provider of genome coordinates
 version           | text                        |           |          |         | extended |              | genome build version for genome coordinates
 _createdby_key    | integer                     |           | not null | 1001    | plain    |              | user who created the record
 _modifiedby_key   | integer                     |           | not null | 1001    | plain    |              | user who last modified the record
 creation_date     | timestamp without time zone |           | not null | now()   | plain    |              | date record was created
 modification_date | timestamp without time zone |           | not null | now()   | plain    |              | date record was last modified
Indexes:
    "mrk_location_cache_pkey" PRIMARY KEY, btree (_marker_key) CLUSTER
    "mrk_location_cache_0" btree (lower(chromosome))
    "mrk_location_cache_idx_chromosome_cmoffset" btree (chromosome, cmoffset)
    "mrk_location_cache_idx_clustered" btree (chromosome, startcoordinate, endcoordinate)
    "mrk_location_cache_idx_marker_type_key" btree (_marker_type_key)
    "mrk_location_cache_idx_organism_key" btree (_organism_key)
Foreign-key constraints:
    "mrk_location_cache__createdby_key_fkey" FOREIGN KEY (_createdby_key) REFERENCES mgi_user(_user_key) DEFERRABLE
    "mrk_location_cache__marker_key_fkey" FOREIGN KEY (_marker_key) REFERENCES mrk_marker(_marker_key) ON DELETE CASCADE DEFERRABLE
    "mrk_location_cache__modifiedby_key_fkey" FOREIGN KEY (_modifiedby_key) REFERENCES mgi_user(_user_key) DEFERRABLE
    "mrk_location_cache__organism_key_fkey" FOREIGN KEY (_organism_key) REFERENCES mgi_organism(_organism_key) DEFERRABLE
*/