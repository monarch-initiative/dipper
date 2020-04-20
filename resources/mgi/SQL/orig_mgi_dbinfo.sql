/*
mgd=> \d+ mgi_dbinfo
                                                           Table "mgd.mgi_dbinfo"
       Column       |            Type             | Collation | Nullable | Default | Storage  | Stats target |          Description          
--------------------+-----------------------------+-----------+----------+---------+----------+--------------+-------------------------------
 public_version     | text                        |           | not null |         | extended |              | public version
 product_name       | text                        |           | not null |         | extended |              | product name
 schema_version     | text                        |           | not null |         | extended |              | schema version
 snp_schema_version | text                        |           | not null |         | extended |              | snp schema version
 snp_data_version   | text                        |           | not null |         | extended |              | snp data version
 lastdump_date      | timestamp without time zone |           | not null | now()   | plain    |              | date of last database dump
 creation_date      | timestamp without time zone |           | not null | now()   | plain    |              | date record was created
 modification_date  | timestamp without time zone |           | not null | now()   | plain    |              | date record was last modified
Indexes:
    "mgi_dbinfo_pkey" PRIMARY KEY, btree (public_version)
*/
