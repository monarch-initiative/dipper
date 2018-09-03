CREATE TABLE `Article_Breed` (
  `article_id` integer  NOT NULL default '0'
,  `breed_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`article_id`,`breed_id`)
,  CONSTRAINT `Article_Breed_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Article_Breed_ibfk_2` FOREIGN KEY (`breed_id`) REFERENCES `Breed` (`breed_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Article_Gene` (
  `article_id` integer  NOT NULL default '0'
,  `gene_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '1'
,  PRIMARY KEY  (`article_id`,`gene_id`)
,  CONSTRAINT `Article_Gene_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Article_Gene_ibfk_4` FOREIGN KEY (`article_id`) REFERENCES `Articles` (`article_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Article_Keyword` (
  `article_id` integer  NOT NULL default '0'
,  `keyword_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`article_id`,`keyword_id`)
,  CONSTRAINT `Article_Keyword_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Article_Keyword_ibfk_2` FOREIGN KEY (`article_id`) REFERENCES `Articles` (`article_id`) ON DELETE CASCADE ON UPDATE CASCADE
);
CREATE TABLE `Article_People` (
  `article_id` integer  NOT NULL default '0'
,  `person_id` integer  NOT NULL default '0'
,  `position` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`article_id`,`person_id`)
,  CONSTRAINT `Article_People_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Article_People_ibfk_2` FOREIGN KEY (`article_id`) REFERENCES `Articles` (`article_id`) ON DELETE CASCADE ON UPDATE CASCADE
);
CREATE TABLE `Article_Phene` (
  `article_id` integer  NOT NULL default '0'
,  `phene_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`article_id`,`phene_id`)
,  CONSTRAINT `Article_Phene_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Article_Phene_ibfk_2` FOREIGN KEY (`phene_id`) REFERENCES `Phene` (`phene_id`) ON DELETE CASCADE ON UPDATE CASCADE
);
CREATE TABLE `Articles` (
  `article_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `title` text NOT NULL
,  `journal` varchar(200) default NULL
,  `volume` varchar(100) default NULL
,  `pages` varchar(105) default NULL
,  `year` year(4) NOT NULL default '0000'
,  `locus` varchar(30) default NULL
,  `abstract` text
,  `publisher` integer  default NULL
,  `pubmed_id` integer  default NULL
,  `library` text
,  `added_by` integer  NOT NULL default '0'
,  `date_modified` date NOT NULL default '0000-00-00'
,  CONSTRAINT `Articles_ibfk_1` FOREIGN KEY (`publisher`) REFERENCES `Publishers` (`publish_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Articles_ibfk_2` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`)
);
CREATE TABLE `Breed` (
  `breed_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `breed_name` varchar(60) NOT NULL default ''
,  `gb_species_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  CONSTRAINT `Breed_ibfk_1` FOREIGN KEY (`gb_species_id`) REFERENCES `Species_gb` (`gb_species_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Breed_ibfk_2` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Breed_Phene` (
  `breed_id` integer  NOT NULL default '0'
,  `phene_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`phene_id`,`breed_id`)
,  CONSTRAINT `Breed_Phene_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Breed_Phene_ibfk_2` FOREIGN KEY (`breed_id`) REFERENCES `Breed` (`breed_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Genes_gb` (
  `gene_id` integer  NOT NULL default '0'
,  `gb_species_id` integer  NOT NULL default '0'
,  `pubmed_id` integer  default '0'
,  `symbol` varchar(45) NOT NULL default ''
,  `gene_desc` text
,  `gene_type` text  NOT NULL default 'protein-coding'
,  `added_by` integer  NOT NULL default '1'
,  `date_modified` date NOT NULL default '0000-00-00'
,  UNIQUE (`gene_id`,`gb_species_id`,`pubmed_id`)
);
CREATE TABLE `Group_Categories` (
  `cat_id` integer  NOT NULL default '0'
,  `cat_name` varchar(100) NOT NULL default ''
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  PRIMARY KEY  (`cat_id`)
,  CONSTRAINT `Group_Categories_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Group_MPO` (
  `omia_id` integer  zerofill NOT NULL default '000000'
,  `MPO_no` integer  zerofill NOT NULL default '000000'
,  `added_by` integer  NOT NULL default '255'
,  PRIMARY KEY  (`omia_id`,`MPO_no`)
,  CONSTRAINT `Group_MPO_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Group_MPO_ibfk_2` FOREIGN KEY (`omia_id`) REFERENCES `OMIA_Group` (`omia_id`) ON DELETE CASCADE ON UPDATE CASCADE
);
CREATE TABLE `Inherit_Type` (
  `inherit_id` char(3) NOT NULL default ''
,  `inherit_name` varchar(50) NOT NULL default ''
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  PRIMARY KEY  (`inherit_id`)
,  UNIQUE (`inherit_name`)
,  CONSTRAINT `Inherit_Type_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Keywords` (
  `keyword_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `keyword` varchar(50) NOT NULL default ''
,  `added_by` integer  NOT NULL default '2'
,  CONSTRAINT `Keywords_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Landmark` (
  `article_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`article_id`)
);
CREATE TABLE `Lida_Links` (
  `lidaurl` varchar(200) NOT NULL default ''
,  `omia_id` integer  zerofill NOT NULL default '000000'
,  `added_by` integer  default NULL
,  PRIMARY KEY  (`omia_id`,`lidaurl`)
,  CONSTRAINT `Lida_Links_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Lida_Links_ibfk_2` FOREIGN KEY (`omia_id`) REFERENCES `OMIA_Group` (`omia_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `OMIA_Group` (
  `omia_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `group_name` varchar(100) NOT NULL default ''
,  `group_summary` text
,  `group_category` integer  default NULL
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  UNIQUE (`omia_id`)
,  CONSTRAINT `OMIA_Group_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `OMIA_Group_ibfk_2` FOREIGN KEY (`group_category`) REFERENCES `Group_Categories` (`cat_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `OMIA_author` (
  `omia_auth_name` varchar(20) NOT NULL default 'data miner'
,  `omia_auth_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  UNIQUE (`omia_auth_id`)
);
CREATE TABLE `Omim_Xref` (
  `omia_id` integer  zerofill NOT NULL default '000000'
,  `omim_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '2'
,  PRIMARY KEY  (`omia_id`,`omim_id`)
,  CONSTRAINT `Omim_Xref_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Omim_Xref_ibfk_2` FOREIGN KEY (`omia_id`) REFERENCES `OMIA_Group` (`omia_id`) ON DELETE CASCADE ON UPDATE CASCADE
);
CREATE TABLE `People` (
  `person_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `name` varchar(100) NOT NULL default ''
,  `added_by` integer  NOT NULL default '2'
,  CONSTRAINT `People_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Phene` (
  `phene_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `omia_id` integer  zerofill NOT NULL default '000000'
,  `gb_species_id` integer  NOT NULL default '0'
,  `phene_name` varchar(60) default NULL
,  `summary` text
,  `symbol` varchar(10) default NULL
,  `marker` text
,  `clin_feat` text
,  `gen_test` text
,  `inherit` char(3) default NULL
,  `inherit_text` text
,  `mol_gen` text
,  `map_info` text
,  `history` text
,  `control` text
,  `pathology` text
,  `prevalence` text
,  `defect` text  NOT NULL default 'no'
,  `singlelocus` text  NOT NULL default 'unknown'
,  `characterised` text  NOT NULL default 'Unknown'
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  CONSTRAINT `Phene_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
,  CONSTRAINT `Phene_ibfk_2` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Phene_Gene` (
  `gene_id` integer  NOT NULL default '0'
,  `phene_id` integer  NOT NULL default '0'
,  `added_by` integer  NOT NULL default '0'
);
CREATE TABLE `Publishers` (
  `publish_id` integer  NOT NULL default '0'
,  `name` varchar(100) default NULL
,  `place` text
,  `added_by` integer  NOT NULL default '2'
,  `date_modified` date NOT NULL default '0000-00-00'
,  PRIMARY KEY  (`publish_id`)
);
CREATE TABLE `Resources` (
  `phene_id` integer  NOT NULL default '0'
,  `resource_url` varchar(60) NOT NULL default ''
,  `resource_id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT
,  `added_by` integer  NOT NULL default '2'
,  CONSTRAINT `Resources_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Species_gb` (
  `gb_species_id` integer  NOT NULL default '0'
,  `sci_name` varchar(60) NOT NULL default ''
,  `com_name` varchar(60) default NULL
,  `added_by` integer  NOT NULL default '1'
,  `date_modified` date NOT NULL default '0000-00-00'
,  UNIQUE (`gb_species_id`)
,  CONSTRAINT `Species_gb_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `Synonyms` (
  `synonym` varchar(100) NOT NULL default ''
,  `accepted_word` varchar(100) NOT NULL default ''
,  `added_by` integer  NOT NULL default '1'
,  `date_modified` date NOT NULL default '0000-00-00'
,  PRIMARY KEY  (`synonym`)
,  CONSTRAINT `Synonyms_ibfk_1` FOREIGN KEY (`added_by`) REFERENCES `OMIA_author` (`omia_auth_id`) ON DELETE NO ACTION ON UPDATE CASCADE
);
CREATE TABLE `tmp` (
  `gene_id` integer  NOT NULL default '0'
,  `gb_species_id` integer  NOT NULL default '0'
,  `pubmed_id` integer  default '0'
,  `symbol` varchar(45) NOT NULL default ''
,  `gene_desc` text
,  `gene_type` text  NOT NULL default 'protein-coding'
,  `added_by` integer  NOT NULL default '1'
);
CREATE INDEX "idx_Omim_Xref_added_by" ON "Omim_Xref" (`added_by`);
CREATE INDEX "idx_Breed_breed_name" ON "Breed" (`breed_name`);
CREATE INDEX "idx_Breed_species_id" ON "Breed" (`gb_species_id`);
CREATE INDEX "idx_Breed_added_by" ON "Breed" (`added_by`);
CREATE INDEX "idx_People_name" ON "People" (`name`);
CREATE INDEX "idx_People_added_by" ON "People" (`added_by`);
CREATE INDEX "idx_Group_MPO_added_by" ON "Group_MPO" (`added_by`);
CREATE INDEX "idx_Article_People_position" ON "Article_People" (`position`);
CREATE INDEX "idx_Article_People_added_by" ON "Article_People" (`added_by`);
CREATE INDEX "idx_Article_People_person_id" ON "Article_People" (`person_id`);
CREATE INDEX "idx_Article_People_article_id" ON "Article_People" (`article_id`);
CREATE INDEX "idx_Phene_phene_id" ON "Phene" (`omia_id`);
CREATE INDEX "idx_Phene_species_id" ON "Phene" (`gb_species_id`);
CREATE INDEX "idx_Phene_symbol" ON "Phene" (`symbol`);
CREATE INDEX "idx_Phene_inherit" ON "Phene" (`inherit`);
CREATE INDEX "idx_Phene_added_by" ON "Phene" (`added_by`);
CREATE INDEX "idx_Synonyms_added_by" ON "Synonyms" (`added_by`);
CREATE INDEX "idx_Phene_Gene_added_by" ON "Phene_Gene" (`added_by`);
CREATE INDEX "idx_Phene_Gene_phene_id" ON "Phene_Gene" (`phene_id`);
CREATE INDEX "idx_Phene_Gene_gene_id" ON "Phene_Gene" (`gene_id`);
CREATE INDEX "idx_Species_gb_added_by" ON "Species_gb" (`added_by`);
CREATE INDEX "idx_Breed_Phene_added_by" ON "Breed_Phene" (`added_by`);
CREATE INDEX "idx_Breed_Phene_breed_id" ON "Breed_Phene" (`breed_id`);
CREATE INDEX "idx_Publishers_name" ON "Publishers" (`name`);
CREATE INDEX "idx_Lida_Links_added_by" ON "Lida_Links" (`added_by`);
CREATE INDEX "idx_Articles_publisher" ON "Articles" (`publisher`);
CREATE INDEX "idx_Articles_year" ON "Articles" (`year`);
CREATE INDEX "idx_Articles_pmid" ON "Articles" (`pubmed_id`);
CREATE INDEX "idx_Articles_added_by" ON "Articles" (`added_by`);
CREATE INDEX "idx_Article_Keyword_added_by" ON "Article_Keyword" (`added_by`);
CREATE INDEX "idx_Article_Keyword_keyword_id" ON "Article_Keyword" (`keyword_id`);
CREATE INDEX "idx_Article_Keyword_article_id" ON "Article_Keyword" (`article_id`);
CREATE INDEX "idx_Group_Categories_cat_name" ON "Group_Categories" (`cat_name`);
CREATE INDEX "idx_Group_Categories_added_by" ON "Group_Categories" (`added_by`);
CREATE INDEX "idx_Article_Breed_added_by" ON "Article_Breed" (`added_by`);
CREATE INDEX "idx_Article_Breed_breed_id" ON "Article_Breed" (`breed_id`);
CREATE INDEX "idx_OMIA_Group_phene_name" ON "OMIA_Group" (`group_name`);
CREATE INDEX "idx_OMIA_Group_added_by" ON "OMIA_Group" (`added_by`);
CREATE INDEX "idx_OMIA_Group_group_category" ON "OMIA_Group" (`group_category`);
CREATE INDEX "idx_Resources_phene_id" ON "Resources" (`phene_id`);
CREATE INDEX "idx_Resources_added_by" ON "Resources" (`added_by`);
CREATE INDEX "idx_Inherit_Type_inherit_name" ON "Inherit_Type" (`inherit_name`);
CREATE INDEX "idx_Inherit_Type_added_by" ON "Inherit_Type" (`added_by`);
CREATE INDEX "idx_Keywords_keyword" ON "Keywords" (`keyword`);
CREATE INDEX "idx_Keywords_added_by" ON "Keywords" (`added_by`);
CREATE INDEX "idx_Landmark_added_by" ON "Landmark" (`added_by`);
CREATE INDEX "idx_Article_Phene_added_by" ON "Article_Phene" (`added_by`);
CREATE INDEX "idx_Article_Phene_phene_id" ON "Article_Phene" (`phene_id`);
CREATE INDEX "idx_Article_Gene_added_by" ON "Article_Gene" (`added_by`);
CREATE INDEX "idx_Article_Gene_gene_id" ON "Article_Gene" (`gene_id`);
/* No STAT tables available */
