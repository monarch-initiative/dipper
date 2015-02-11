import csv
import os
from datetime import datetime
from stat import *
import re
import psycopg2


from rdflib import Literal
from rdflib.namespace import RDFS, OWL, RDF, DC
from rdflib import Namespace, URIRef, BNode, Graph

from utils import pysed
from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from utils.CurieUtil import CurieUtil
import config
import curie_map
from utils.GraphUtils import GraphUtils


class MGI(Source):
    '''
    Be sure to have pg user/password connection details in your conf.json file, like:
      dbauth : {
        'mgi' : {'user' : '<username>', 'password' : '<password>'}
      }
    '''
# tables in existing interop

    #USED
    # gxd_genotype_view, gxd_genotype_summary_view, gxd_allelepair_view, all_summary_view,
    # all_allele_view, all_allele_mutation_view, mrk_marker_view, voc_annot_view,
    # voc_evidence_view, bib_acc_view, prb_strain_view

    #NOT YET USED
    # mgi_organism_acc_view: Don't think we need this as I have handled the taxon mapping through a map_taxon function, unless we want the MGI ID for the organism.
    # mgi_organism_view: I mapped the taxon from the mrk_marker_view to include all used organisms, but the mapping could be done in a more complete fashion with this table.
    # mgi_reference_allele_view: Don't believe this view is used in either the genotype of phenotype view
    # all_allele_cellline_view: Don't believe this view is used in either the genotype of phenotype view
    # voc_term_view: Don't believe this view is used in either the genotype of phenotype view
    # mrk_summary_view: Used in genotype view. Only need it if we want the MGI ID for the gene.
    # mgi_note_vocevidence_view: Used in phenotype view for free_text_phenotype_description
    # acc_logicaldb_view: Don't believe this view is used in either the genotype of phenotype view
    # mgi_note_strain_view: Don't believe this view is used in either the genotype of phenotype view
    # FIXME - prb_strain_acc_view: Don't believe this view is used in either the genotype of phenotype view, but is needed if we want the MGI ID for the strain.
    # prb_strain_summary_view: Don't believe this view is used in either the genotype of phenotype view
    # prb_strain_marker_view: Don't believe this view is used in either the genotype of phenotype view

#TODO: QA List
    #Do identifiers have the proper prefixes?
    #Are there quotes that need to be stripped from variables?
    #Is there scrubbing needed for any variables?
    #If present in other functions, can the scrubbing be moved to the scrub function?
    #Make a checklist for the full graph and confirm that all nodes are present.
    #Do we need to do any HTML formatting of labels? (< -> &lt;)

    tables = [
        'mgi_dbinfo',
        'gxd_genotype_view',
        'gxd_genotype_summary_view',
        'gxd_allelepair_view',
        'all_summary_view',
        'all_allele_view',
        'all_allele_mutation_view',
        'mrk_marker_view',
        'voc_annot_view',
        'voc_evidence_view',
        'bib_acc_view',
        'prb_strain_view',
        'mrk_summary_view',
        'mrk_acc_view',
        'prb_strain_acc_view'
    ]


    relationship = {
        'is_mutant_of' : 'GENO:0000440',
        'derives_from' : 'RO:0001000',
        'has_alternate_part' : 'GENO:0000382',
        'has_reference_part' : 'GENO:0000385',
        'in_taxon' : 'RO:0000216',
        'has_zygosity' : 'GENO:0000608',
        'is_sequence_variant_instance_of' : 'GENO:0000408',
        'is_reference_instance_of' : 'GENO:0000610',
        'hasExactSynonym' : 'OIO:hasExactSynonym',
        'has_disposition' : 'GENO:0000208',
        'has_phenotype' : 'RO:0002200',
        'has_part' : 'BFO:0000051'
    }

    terms = {
        'variant_locus' : 'GENO:0000002',
        'reference_locus' : 'GENO:0000036',
        'sequence_alteration' : 'SO:0001059',
        'variant_single_locus_complement' : 'GENO:0000030',
        'allele' : 'GENO:0000008',
        'intrinsic_genotype' : 'GENO:0000000',
        'phenotype' : 'MONARCH:phenotype',  # Is this correct? What about GENO:0000348 - phenotype? MONARCH:phenotype
        'evidence' : 'MONARCH:evidence',
        'genomic_background' : 'GENO:0000010'
    }

    def __init__(self):
        Source.__init__(self, 'mgi')
        self.namespaces.update(curie_map.get())
        #assemble all the curie mappings from the imported models
        #self.namespaces.update(Assoc.curie_map)
        #self.namespaces.update(Genotype.curie_map)
        #self.namespaces.update(G2PAssoc.curie_map)

        #update the dataset object with details about this resource
        self.dataset = Dataset('mgi', 'MGI', 'http://www.informatics.jax.org/')

        #check if config exists; if it doesn't, error out and let user know
        if (not (('dbauth' in config.get_config()) and ('mgi' in config.get_config()['dbauth']))):
            print("ERROR: not configured with PG user/password.")

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return

    def fetch(self, is_dl_forced):
        '''
        For the MGI resource, we connect to the remote database, and pull the tables into local files.
        We'll check the local table versions against the remote version
        :return:
        '''

        #create the connection details for MGI
        cxn = config.get_config()['dbauth']['mgi']
        cxn.update({'host' : 'adhoc.informatics.jax.org', 'database' : 'mgd', 'port' : 5432 })

        self.dataset.setFileAccessUrl(('').join(('jdbc:postgresql://',cxn['host'],':',str(cxn['port']),'/',cxn['database'])))

        #process the tables
        #self.fetch_from_pgdb(self.tables,cxn,100)  #for testing
        self.fetch_from_pgdb(self.tables,cxn)

        datestamp=ver=None
        #get the resource version information from table mgi_dbinfo, already fetched above
        outfile=('/').join((self.rawdir,'mgi_dbinfo'))

        if os.path.exists(outfile):
            st = os.stat(outfile)
            with open(outfile, 'r') as f:
                f.readline() #read the header row; skip
                info = f.readline()
                cols = info.split('\t')
                ver = cols[0] #col 0 is public_version
                ver = ver.replace('MGI ','')  #MGI 5.20 --> 5.20
                #MGI has a datestamp for the data within the database; use it instead of the download date
                #datestamp in the table: 2014-12-23 00:14:20
                d = cols[7].strip()  #modification date
                datestamp = datetime.strptime(d, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%d")
                f.close()

        self.dataset.setVersion(datestamp,ver)

        return

    def scrub(self):
        #TODO any scrubbing needed for this resource?
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes: (none)
        :return: None
        '''

        #Should the wild type alleles be removed?

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")

        #so that we don't have to deal with BNodes, we will create hash lookups for the internal identifiers
        #the hash will hold the type-specific-object-keys to MGI public identifiers.  then, subsequent
        #views of the table will lookup the identifiers in the hash.  this allows us to do the 'joining' on the
        #fly

        self.idhash = {'allele' : {}, 'marker' : {}, 'publication' : {}, 'strain' : {}, 'genotype' : {}, 'annot' : {}}

        #the following will provide us the hash-lookups
        self._process_mrk_acc_view(('/').join((self.rawdir,'mrk_acc_view')),limit) #DONE
        self._process_all_summary_view(('/').join((self.rawdir,'all_summary_view')),limit) #DONE
        self._process_bib_acc_view(('/').join((self.rawdir,'bib_acc_view')),limit)  #DONE
        self._process_gxd_genotype_summary_view(('/').join((self.rawdir,'gxd_genotype_summary_view')),limit)  #DONE
        self._process_prb_strain_acc_view(('/').join((self.rawdir,'prb_strain_acc_view')),limit)  #DONE
        #FIXME add prb_strain_acc_view

        #the following will use the hash to lookup the ids
        self._process_prb_strain_view(('/').join((self.rawdir,'prb_strain_view')),limit)
        self._process_gxd_genotype_view(('/').join((self.rawdir,'gxd_genotype_view')),limit)  #DONE
        self._process_all_allele_view(('/').join((self.rawdir,'all_allele_view')),limit) #DONE
        self._process_gxd_allele_pair_view(('/').join((self.rawdir,'gxd_allelepair_view')),limit)  #DONE
        self._process_all_allele_mutation_view(('/').join((self.rawdir,'all_allele_mutation_view')),limit) #DONE
        self._process_mrk_summary_view(('/').join((self.rawdir,'mrk_summary_view')),limit)  #DONE
        self._process_mrk_marker_view(('/').join((self.rawdir,'mrk_marker_view')),limit)  #DONE
        self._process_voc_annot_view(('/').join((self.rawdir,'voc_annot_view')),limit)  #DONE
        self._process_voc_evidence_view(('/').join((self.rawdir,'voc_evidence_view')),limit)  #DONE


        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Loaded", len(self.graph), "nodes")
        return


    def _process_gxd_genotype_view(self, raw, limit=None):
        '''
        This table indicates the relationship between a genotype and it's background strain.

        Makes these triples:
        <MGI:genotypeid> GENO:has_reference_part <internal strain id>

        if the genotype id isn't in the hashmap, it adds it here (but this shouldn't happen):
        <MGI:genotypeid> a GENO:genotype
        <MGI:genotypeid> sameAs <internal genotype id>
        <internal strain id> a GENO:genomic_background

        :param raw:
        :param limit:
        :return:
        '''

        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f1:
            f1.readline()  # read the header row; skip
            for line in f1:
                line_counter += 1
                (genotype_key,strain_key,isconditional,note,existsas_key,createdby_key,modifiedby_key,creation_date,
                 modification_date,strain,mgiid,dbname,createdbymodifiedby,existsas,empty) = line.split('\t')

                if self.idhash['genotype'].get(genotype_key) is None:
                    #just in case we haven't seen it before, catch and add the id mapping here
                    self.idhash['genotype'][genotype_key] = mgiid
                    gu.addIndividualToGraph(self.graph,mgiid,None,self.terms['intrinsic_genotype'])
                    #TODO get label
                gt = URIRef(cu.get_uri(mgiid))

                #if it's in the hash, assume that the individual was created elsewhere
                strain_id = self.idhash['strain'].get(strain_key)
                if (strain_id is not None):
                    str = URIRef(cu.get_uri(strain_id))

                   #todo refactor to genotype model
                    self.graph.add((gt,URIRef(cu.get_uri(self.relationship['has_reference_part'])),str))

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_gxd_genotype_summary_view(self,raw,limit=None):
        #need to make triples:
        #. genotype is a class - redundant?
        #. genotype has equivalent class internalGenotypeID
        #. genotype subclass of intrinsic_genotype - redundant?
        #. genotype has label description

        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,createdby_key,modifiedby_key,
                 creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')

                #add the internal genotype to mgi mapping
                self.idhash['genotype'][object_key] = mgiid

                #FIXME note the short_description is the GVC  (use this or reason?)

                if (preferred == '1'):
                    gu.addIndividualToGraph(self.graph,mgiid,description,self.terms['intrinsic_genotype'])
                #TODO what to do with != preferred

                if (limit is not None and line_counter > limit):
                    break

        return

    #NOTE: might be best to process alleles initially from the all_allele_view, as this does not have any repeats of alleles!
    def _process_all_summary_view(self,raw,limit):
        '''
        Here, we get the allele definitions: id, label, description, type
        Also, we add the id to this source's global hash for lookup later
        :param raw:
        :param limit:
        :return:
        '''
        #Need to make triples:
        #. allele is an instance of allele
        #. internalAlleleID has equivalent class as allele
        #. allele has label short_description: Better to use symbol from all_allele_view?
        #. allele has description description

        #TODO: allele subtype

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                 createdby_key,modifiedby_key,creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')
                #NOTE:May want to filter alleles based on the preferred field (preferred = 1) or will get duplicates
                ## (24288, to be exact... Reduced to 480 if filtered on preferred = 1)

                #If we want to filter on preferred:
                if preferred == '1':
                    #add the allele key to the hash for later lookup
                    self.idhash['allele'][object_key] = mgiid
                    gu.addIndividualToGraph(self.graph,mgiid,short_description.strip(),self.terms['allele'],description.strip())
                    #TODO update subtype

                #TODO deal with non-preferreds, are these deprecated?

                if (limit is not None and line_counter > limit):
                    break

        return


    def _process_all_allele_view(self,raw,limit):
        """
        Add the allele as a variant locus (or reference locus if wild-type).
        If the marker is specified, we add the link to the marker.
        We assume that the MGI ids are available in the idhash, added in all_summary_view.
        We add the sequence alteration as a BNode here.

        Triples:
        <MGI:allele_id> a GENO:variant_locus OR GENO:reference_locus
        <MGI:allele_id> GENO:has_variant_part OR GENO:has_reference_part  <MGI:marker_id>
        <MGI:allele_id> GENO:derived_from <MGI:strain_id>
        <MGI:allele_id>
        :param raw:
        :param limit:
        :return:
        """

        #NOTE: allele == variant locus
        #Need triples:
        #. (variant) allele is a subclass of variant_locus
        #. (variant) allele is variant_of gene/marker
        #. (wild type) allele is a subclass of reference_locus
        #. (wild type) allele is reference_of gene/marker
        #. allele has label symbol (any reformatting?)
        #. sequence alteration is a class
        #. sequence alteration is a subclass of SO:0001059
        #. sequence alteration has description name
        #. sequence alteration in strain

        # Extra: strain_key, map along the lines of "allele (allele_key -> Bnode) in strain (strain_key -> Bnode)?"
        # Strain label available. Marker label available. Better to map those through their primary tables, correct?
        #TODO
        # Allele type key also available. Need to locate related table
        # transmission_key -> inheritance? Need to locate related table.
        # strain: sequence_alteration in strain?

        #Instead of a function-specific set of variables, should these instead be added
        # to the relationship table at the top?

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        print("INFO: adding alleles; mapping to markers; extracting their sequence alterations")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allele_key,marker_key,strain_key,mode_key,allele_type_key,allele_status_key,transmission_key,
                 collection_key,symbol,name,nomensymbol,iswildtype,isextinct,ismixed,createdby_key,modifiedby_key,
                 approvedby_key,approval_date,creation_date,modification_date,markersymbol,term,statusnum,strain,collection,createdby,modifiedby,approvedby) = line.split('\t')

                allele_id = self.idhash['allele'].get(allele_key)
                if (allele_id is None):
                    print("ERROR: what to do! can't find allele_id. skipping",allele_key,symbol)
                    continue

                if (marker_key is not None) and (marker_key != ''):
                    #we make the assumption here that the markers have already been added to the table
                    marker_id = self.idhash['marker'].get(marker_key)
                    if (marker_id is None):
                        print("ERROR: what to do! can't find marker_id. skipping",marker_key,symbol)
                        continue

                #TODO handle the case when the marker == allele (though they are not the same MGI id)

                strain_id = self.idhash['strain'].get(strain_key)
                iseqalt_id = self._makeInternalIdentifier('seqalt',allele_key)
                iseqalt = BNode(iseqalt_id)  # Any issues with reusing the allele_key as long as we use a different prefix?

                # for non-wild type alleles:
                if iswildtype == '0':
                    locus_type = self.terms['variant_locus']
                    rel = URIRef(cu.get_uri(self.relationship['is_sequence_variant_instance_of']))
                #for wild type alleles:
                elif iswildtype == '1':
                    locus_type = self.terms['reference_locus']
                    rel = URIRef(cu.get_uri(self.relationship['is_reference_instance_of']))

                gu.addIndividualToGraph(self.graph,allele_id,symbol,locus_type)
                #add link between gene and allele
                al = URIRef(cu.get_uri(allele_id))
                if marker_id is not None:
                    #marker_id will be none if the allele is not linked to a marker (as in, it's not mapped to a locus)
                    self.graph.add((al,
                                    rel,
                                    URIRef(cu.get_uri(marker_id))))



                #sequence alteration in strain
                #FIXME change this to a different relation in_strain, genomically_related_to, sequence_derives_from
                if iswildtype == '0':
                    #sequence alteration has label reformatted(symbol)
                    sa_label = symbol
                    if re.match(".*<.*>.*", symbol):
                        #print(sa_label)
                        sa_label = re.sub(".*<", "<", symbol)
                        #print(sa_label)
                    elif re.match("\+", symbol):
                        #TODO: Check to see if this is the proper handling, as while symbol is just +, marker symbol has entries without any <+>.
                        sa_label = '<+>'
                        #print(sa_label)
                    #removing the < and > from sa
                    sa_label = re.sub('[\<\>]','',sa_label)
                    gu.addIndividualToGraph(self.graph,iseqalt_id,sa_label,self.terms['sequence_alteration'],name)
                    rel = URIRef(cu.get_uri(self.relationship['has_alternate_part']))
                    self.graph.add((al,rel,iseqalt))  #TODO IS THIS RIGHT?


                    if strain_id is not None:
                        self.graph.add((al,URIRef(cu.get_uri(self.relationship['derives_from'])),URIRef(cu.get_uri(strain_id))))

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_gxd_allele_pair_view(self,raw,limit):
        #Need triples:
        #. vslc is of type: vslc
        #. genotype has vslc allele_pair_key
        #. vslc has label processed(vslc_label)
        #. vslc has_part allele1
        #. vslc has_part allele2
        #. vslc has_disposition mapped(allelestate)


        #Additional stuff: chromosome, compound? (entries: Top, Not Applicable, Bottom)

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allelepair_key,genotype_key,allele_key_1,allele_key_2,marker_key,mutantcellline_key_1,mutantcellline_key_2,
                 pairstate_key,compound_key,sequencenum,createdby_key,modifiedby_key,creation_date,modification_date,symbol,
                 chromosome,allele1,allele2,allelestate,compound) = line.split('\t')
                #NOTE: symbol = gene/marker, allele1 + allele2 = VSLC, allele1/allele2 = variant locus, allelestate = zygosity
                #FIXME Need to handle alleles not in the *<*> format, such as many gene traps, induced mutations, and transgenics

                genotype_id = self.idhash['genotype'].get(genotype_key)

                allele1_id = self.idhash['allele'].get(allele_key_1)
                allele2_id = self.idhash['allele'].get(allele_key_2)

                #Need to map the allelestate to a zygosity term
                zygosity_id = self._map_zygosity(allelestate)
                ivslc = BNode('vslckey'+allelepair_key)
                #FIXME: VSLC label likely needs processing similar to the processing in the all_allele_view
                #FIXME: Handle null alleles for allele2
                vslc_label = (allele1+'/'+allele2)
                #print(vslc_label)

                #. vslc is of type: vslc
                self.graph.add((ivslc,RDF['type'],URIRef(cu.get_uri(self.terms['variant_single_locus_complement']))))

                #. vslc has label processed(vslc_label)
                self.graph.add((ivslc,RDFS['label'],Literal(vslc_label)))

                #genotype has part vslc
                self.graph.add((URIRef(cu.get_uri(genotype_id)),URIRef(cu.get_uri(self.relationship['has_alternate_part'])),ivslc))

                #vslc has parts allele1/allele2
                rel = URIRef(cu.get_uri(self.relationship['has_part']))
                if allele1_id is not None:
                    self.graph.add((ivslc,rel,URIRef(cu.get_uri(allele1_id))))
                if allele2_id is not None:
                    self.graph.add((ivslc,rel,URIRef(cu.get_uri(allele2_id))))

                #FIXME: Is this correct?
                # Also, in my concept map I had zygosity as GENO:0000608 - has_zygosity,
                # but I don't see it in my geno.owl file.
                self.graph.add((ivslc,URIRef(cu.get_uri(self.relationship['has_zygosity'])),URIRef(cu.get_uri(zygosity_id))))

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_all_allele_mutation_view(self,raw,limit):
        #Need triples:
        #. sequence_alteration has_type mutation
        #. sequence_alteration_type_label?

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (allele_key,mutation_key,creation_date,modification_date,mutation) = line.split('\t')
                iseqalt_id = self._makeInternalIdentifier('seqalt',allele_key)
                iseqalt = BNode(iseqalt_id)

                #map the sequence_alteration_type
                seq_alt_type_id = self._map_seq_alt_type(mutation)
                self.graph.add((iseqalt,RDF['type'],URIRef(cu.get_uri(seq_alt_type_id))))

                if (limit is not None and line_counter > limit):
                    break

        return



    def _process_voc_annot_view(self,raw,limit):
        '''
        This MGI table represents associations between things.
        We currently filter this table on Genotype-Phenotype associations, but may be expanded in the future.
        Reference the genotype in the idhash.
        Add the internal annotation id to the idhashmap

        :param raw:
        :param limit:
        :return:
        '''

        #TODO also get Strain/Attributes (annottypekey = 1000)
        #TODO what is Phenotype (Derived) vs non-derived?  (annottypekey = 1015)
        #TODO is evidence in this table?  what is the evidence vocab key?

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:

                (annot_key,annot_type_key,object_key,term_key,qualifier_key,creation_date,modification_date,qualifier,
                 term,sequence_num,accid,logicaldb_key,vocab_key,mgi_type_key,evidence_vocab_key,anot_type) = line.split('\t')


                # Restricting to type 1002, as done in the MousePhenotypes view.
                # Corresponds to 'Mammalian Phenotype/Genotype' and MP terms
                if annot_type_key == '1002':
                    line_counter += 1

                    #todo add NOT annotations
                    #skip 'normal'
                    if (qualifier=='norm'):
                        print("INFO: found normal phenotype:",term)
                        continue

                    #. This is the phenotype, or MP term for the phenotype.
                    gu.addClassToGraph(self.graph,accid,None)

                    iassoc_id = self._makeInternalIdentifier('annot',annot_key)
                    assoc_id = ':'+self.make_id(iassoc_id)
                    #add the assoc to the hashmap (using the monarch id)
                    self.idhash['annot'][annot_key] = assoc_id
                    genotype_id = self.idhash['genotype'].get(object_key)
                    if (genotype_id is None):
                        print("ERROR: can't find genotype id from",object_key)
                    else:
                        #add the association
                        assoc = G2PAssoc(assoc_id,genotype_id,accid,None,None)
                        assoc.addAssociationNodeToGraph(self.graph)


                if (limit is not None and line_counter > limit):
                    break

        return


    def _process_voc_evidence_view(self,raw,limit):
        """
        Here we fetch the evidence (code and publication) for the associations.
        We will only add the evidence if the annotation is in our idhash.
        :param raw:
        :param limit:
        :return:
        """

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (annot_evidence_key,annot_key,evidence_term_key,refs_key,inferred_from,created_by_key,modified_by_key,
                creation_date,modification_date,evidence_code,evidence_seq_num,jnumid,jnum,short_citation,created_by,modified_by)= line.split('\t')

                assoc_id = self.idhash['annot'].get(annot_key)

                if (assoc_id is None):
                    #assume that we only want to add the evidence/source for annots that we have in our db
                    continue

                # Only 18 evidence codes used in MGI, so create a mapping function to map the label and the ID.
                evidence_id = self._map_evidence_id(evidence_code)
                evidence = URIRef(cu.get_uri(evidence_id))

                #TODO add it as an instance of what type?
                #add the pub as an individual;
                gu.addIndividualToGraph(self.graph,jnumid,None)
                pub = URIRef(cu.get_uri(jnumid))

                #add the ECO and citation information to the annot
                self.graph.add((URIRef(cu.get_uri(assoc_id)),DC['evidence'],evidence))
                self.graph.add((URIRef(cu.get_uri(assoc_id)),DC['source'],pub))


                if (limit is not None and line_counter > limit):
                    break

        return


    def _process_bib_acc_view(self,raw,limit):
        '''
        This traverses the table twice; once to look up the internal key to J number mapping; then to make the
        equivalences.  All internal keys have a J and MGI identifier.
        This will make equivalences between the different pub ids
        :param raw:
        :param limit:
        :return:
        '''



        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)

        #firstpass, get the J number mapping, and add to the global hash
        line_counter = 0
        print('INFO: getting pub id mappings')
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue #skip header
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                created_by_key,modified_by_key,creation_date,modification_date,logical_db)= line

                #we use the J number here because it is the externally-accessible identifier
                if prefixpart != 'J:':
                    continue
                self.idhash['publication'][object_key] = accid
                gu.addIndividualToGraph(self.graph,accid,None)

                if (limit is not None and line_counter > limit):
                    break


        #2nd pass, look up the MGI identifier in the hash
        print("INFO: getting pub equivalent ids")
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1
                if (line_counter == 1):
                    continue #skip header
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,preferred,
                created_by_key,modified_by_key,creation_date,modification_date,logical_db)= line


                logical_db = logical_db.strip()

                #get the real nice pub identifier
                jid = self.idhash['publication'].get(object_key)

                pub_id = None
                if (logicaldb_key == '29'):  #pubmed
                    pub_id = 'PMID:'+accid
                elif (logicaldb_key == '1' and re.match('MGI:',prefixpart)):
                    #don't get the J numbers, because we dont' need to make the equiv to itself.
                    pub_id = accid
                elif (logical_db == 'Journal Link'):
                    #some DOIs seem to have spaces
                    #FIXME MGI needs to FIX THESE UPSTREAM!!!!
                    #we'll scrub them here for the time being
                    pub_id = 'DOI:'+re.sub('\s+','',accid)
                elif (logicaldb_key == '1' and re.match('J:',prefixpart)):
                    #we can skip the J numbers
                    continue

                if (pub_id is not None):
                    #only add these to the graph if it's mapped to something we understand
                    gu.addIndividualToGraph(self.graph,pub_id,None)
                    #todo add this to graph utils
                    self.graph.add((URIRef(cu.get_uri(jid)),OWL['sameAs'],URIRef(cu.get_uri(pub_id))))

                else:
                    print("WARN: Publication from (", logical_db, ") not mapped for",object_key)

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_prb_strain_view(self,raw,limit):
        '''
        Process a table to get strains (with internal ids), and their labels.
        These strains are created as instances of intrinsic_genotype.

        :param raw:
        :param limit:
        :return:
        '''
        #Only 9 strain types if we want to map them (recombinant congenci, inbred strain, NA, congenic,
        # consomic, coisogenic, recombinant inbred, NS, conplastic)
        #160 species types, but could probably slim that down.
        #If we don't want anything else from this table other than the strain label, could potentially drop it
        # and just keep the strain labelling in the gxd_genotype_view.


        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0

        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for line in filereader:
                line_counter += 1

                (strain_key,species_key,strain_type_key,strain,standard,private,genetic_background,created_by_key,
                modified_by_key,creation_date,modification_date,species,strain_type,created_by,modified_by) = line

                strain_id = self.idhash['strain'].get(strain_key)

                istrain_id = self._makeInternalIdentifier('strain',strain_key)

                #FIXME is the strain an 'intrinsic_genotype', 'genomic_background' or something else?
                gu.addIndividualToGraph(self.graph,istrain_id,strain,self.terms['intrinsic_genotype'])
                #TODO add species
                #TODO what is strain type anyway?
                #ispecies = BNode('specieskey'+species_key)
                #istrain_type = BNode('straintypekey'+strain_type_key)

                if (limit is not None and line_counter > limit):
                    break

        return


    def _process_mrk_marker_view(self,raw,limit):
        '''
        This is the definition of markers (as in genes, but other genomic loci as well).
        It looks up the identifiers in the hashmap
        This includes their labels, specific class, and identifiers
        TODO should we use the mrk_mouse_view instead?
        :param raw:
        :param limit:
        :return:
        '''
        #Need triples:
        #. marker is type class
        #. marker has subclass mapped(markertype)
        #. marker has label symbol
        #. marker has synonym name
        #. or marker has description name?
        #Process based on status? (official, withdrawn, interim)
        #Do we want the chromosome number?

        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                (marker_key,organism_key,marker_status_key,marker_type_key,curationstate_key,symbol,name,chromosome,
                cytogenetic_offset,createdby_key,modifiedby_key,creation_date,modification_date,organism,common_name,
                latin_name,status,marker_type,curation_state,created_by,modified_by) = line.split('\t')

                #Remove the withdrawn markers
                if marker_status_key != '2':
                    marker_id = self.idhash['marker'].get(marker_key)
                    #only pull info for mouse genes here?  other species should come from other dbs
                    if (organism_key != '1'):
                        continue

                    if (marker_id is None):
                        print("ERROR: can't find",marker_key,symbol,"in the id hash")

                    #map the marker to the gene class
                    mapped_marker_type = self._map_marker_type(marker_type)

                    gu.addClassToGraph(self.graph,marker_id,symbol,mapped_marker_type,name)
                    gu.addSynonym(self.graph,marker_id,name,Assoc.relationships['hasExactSynonym'])

                    #add the taxon
                    taxon_id = self._map_taxon(latin_name)
                    self.graph.add((URIRef(cu.get_uri(marker_id)),URIRef(cu.get_uri(self.relationship['in_taxon'])),URIRef(cu.get_uri(taxon_id))))

                    #TODO: Think it would make more sense to map the taxon using one of the organism tables.
                    #TODO: If mapping to taxon using an organism table, map to the organism BNode
                    #self.graph.add((imarker,URIRef(cu.get_uri(self.relationship['in_taxon'])),iorganism))

                    if (limit is not None and line_counter > limit):
                        break

        return


    def _process_mrk_summary_view(self,raw,limit):
        #NOTE: There are multiple identifiers available for markers/genes from 28 different resources in this table.
        #Currently handling identifiers from TrEMBL, PDB, ENSEMBL, PRO, miRBASE, MGI, Entrez gene, RefSeq,
        # swiss-prot, and EC, but can add more.

        #Need to grab the iMarker ID, MGI ID
        #Determine if the row is the MGI ID row
            #Process differently if it is. Add to graph as URI?
        #Otherwise process as a non MGI ID row
        #Is it from one of the resources that you wish to use?
            #If so, add as imarker same as marker ID?

        #Need triples:
        #.


        gu = GraphUtils(self.namespaces)
        cu = CurieUtil(self.namespaces)
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1

                #print(line.split('\t'))
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,mgiid,subtype,description,short_description) = line.split('\t')

                if self.idhash['marker'].get(object_key) is None:
                    #can't find the marker in the hash; add it here:
                    self.idhash['marker'][object_key] = mgiid
                    gu.addClassToGraph(self.graph,mgiid,None)

                marker = URIRef(cu.get_uri(mgiid))

                if (accid == mgiid):
                    #don't need to make equivlances to itself
                    continue

                #ISSUE: MirBase accession ID can map to multiple MGI IDs if the miRNA is also part of a cluster (Mirlet7b is part of cluster Mirc31)

                if (preferred=='1'):
                    #logicaldb_key in ['41', '45', '60', '133', '134', '135', '27', '83', '1', '55', '13', '8']: # '27'

                    #Do something
                    mapped_id = None
                    if logicaldb_key in ['133','134','60']:
                        #FIXME some ensembl ids are genes, others are proteins
                        mapped_id = 'ENSEMBL:'+accid
                    elif logicaldb_key == '83':
                        mapped_id = 'miRBase:'+accid
                    elif logicaldb_key == '1':
                        mapped_id = 'MGI:'+accid
                    elif logicaldb_key == '41':
                        mapped_id = 'TrEMBL:'+accid
                    elif logicaldb_key == '45':
                        mapped_id = 'PDB:'+accid
                    elif logicaldb_key == '135':
                        mapped_id = 'PR:'+accid
                    elif logicaldb_key == '83':
                        mapped_id = 'miRBase:'+accid
                    elif logicaldb_key == '55':
                        mapped_id = 'NCBIGene:'+accid
                    elif logicaldb_key == '27':
                        mapped_id = 'RefSeq:'+accid
                    elif logicaldb_key == '13':
                        mapped_id = 'SwissProt:'+accid
                    elif logicaldb_key == '8':
                        mapped_id = 'EC:'+accid
                    #FIXME: The EC IDs are used for multiple genes, resulting in one EC number

                    if (mapped_id is not None):
                        gu.addClassToGraph(self.graph,mapped_id,None)
                        gu.addEquivalentClass(self.graph,mgiid,mapped_id)


                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_mrk_acc_view(self,raw,limit):
        '''
        Use this table to create the idmap between the internal marker id and the public mgiid.
        Also, add the equivalence statements between genes for a subset of the identifiers
        (thus far is NCBIGene and ENSEMBL)
        :param raw:
        :param limit:
        :return:
        '''

        #make a pass through the table first, to create the mapping between the external and internal identifiers
        line_counter = 0
        gu = GraphUtils(self.namespaces)
        print("INFO: mapping markers to internal identifiers")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefix_part,numeric_part,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,logicaldb,organism_key) = line.split('\t')

                #get the hashmap of the identifiers
                if (logicaldb_key == '1') and (prefix_part == 'MGI:') and (preferred == '1'):
                    self.idhash['marker'][object_key] = accid
                    gu.addClassToGraph(self.graph,accid,None)

        #pass through the file again, and make the equivalence statements to a subset of the idspaces
        print("INFO: mapping marker equivalent identifiers")
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefix_part,numeric_part,logicaldb_key,object_key,mgi_type_key,private,preferred,
                 created_by_key,modified_by_key,creation_date,modification_date,logicaldb,organism_key) = line.split('\t')


                mgiid = self.idhash['marker'].get(object_key)
                if (mgiid is None):
                    #presumably we've already added the relevant MGI ids already, so skip those that we can't find
                    #print("INFO:can't find mgiid for",object_key)
                    continue
                marker_id = None
                if (preferred == '1'):  #what does it mean if it's 0?
                    if logicaldb_key == '55':  #entrez/ncbi
                        marker_id = 'NCBIGene:'+accid
                    elif logicaldb_key == '1' and prefix_part != 'MGI:':  #mgi
                        marker_id = accid
                    elif logicaldb_key == '60':
                        marker_id = 'ENSEMBL:'+accid
                    #TODO get other identifiers
                    #TODO get non-preferred ids==deprecated?

                if (marker_id is not None):
                    gu.addClassToGraph(self.graph,marker_id,None)
                    gu.addEquivalentClass(self.graph,mgiid,marker_id)

                if (limit is not None and line_counter > limit):
                    break

        return

    def _process_prb_strain_acc_view(self,raw,limit):
        '''
        Use this table to create the idmap between the internal marker id and the public mgiid.
        Also, add the equivalence statements between genes for a subset of the identifiers
        (thus far is NCBIGene and ENSEMBL)
        :param raw:
        :param limit:
        :return:
        '''

        #make a pass through the table first, to create the mapping between the external and internal identifiers
        line_counter = 0
        gu = GraphUtils(self.namespaces)
        print("INFO: mapping markers to internal identifiers")
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,
                 preferred,createdby_key,modifiedby_key,creation_date,modification_date,logicaldb) = line.split('\t')

                #get the hashmap of the identifiers
                if (logicaldb_key == '1') and (prefixpart == 'MGI:') and (preferred == '1'):
                    self.idhash['strain'][object_key] = accid
                    gu.addClassToGraph(self.graph,accid,None)

        #pass through the file again, and make the equivalence statements to a subset of the idspaces
        print("INFO: mapping strain equivalent identifiers")
        line_counter = 0
        with open(raw, 'r') as f:
            f.readline()  # read the header row; skip
            for line in f:
                line_counter += 1
                (accession_key,accid,prefixpart,numericpart,logicaldb_key,object_key,mgitype_key,private,
                 preferred,createdby_key,modifiedby_key,creation_date,modification_date,logicaldb) = line.split('\t')


                mgiid = self.idhash['strain'].get(object_key)
                if (mgiid is None):
                    #presumably we've already added the relevant MGI ids already, so skip those that we can't find
                    #print("INFO:can't find mgiid for",object_key)
                    continue
                strain_id = None
                if (preferred == '1'):  #what does it mean if it's 0?
                    if logicaldb_key == '22':  #JAX
                        #scrub out the backticks from accids
                        #TODO notify the source upstream
                        accid = re.sub('`','',accid)
                        strain_id = 'JAX:'+accid
                    #TODO get non-preferred ids==deprecated?

                if (strain_id is not None):
                    gu.addClassToGraph(self.graph,strain_id,None)
                    gu.addEquivalentClass(self.graph,mgiid,strain_id)

                if (limit is not None and line_counter > limit):
                    break

        return


    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


    #TODO generalize this to a set of utils
    def _getcols(self,cur,table):
        query=(' ').join(("SELECT * FROM",table,"LIMIT 0"))  #for testing
        #print("COMMAND:",query)
        cur.execute(query)
        colnames = [desc[0] for desc in cur.description]
        print("COLS ("+table+"):",colnames)

        return


    def file_len(self,fname):
        with open(fname) as f:
            return sum(1 for line in f)


    #TODO: Finish identifying SO/GENO terms for mappings for those found in MGI
    def _map_seq_alt_type(self, sequence_alteration_type):
        type = None
        type_map = {
            'Deletion': 'SO:0000159',  # deletion
            'Disruption caused by insertion of vector': 'SO:0000667',  # insertion - correct?
            'Duplication': 'SO:1000035',  # duplication
            'Insertion': 'SO:0000667',  # insertion
            'Insertion of gene trap vector': 'SO:0000667',  # insertion - correct?
            'Intergenic deletion': 'SO:0000159',  # deletion
            'Intragenic deletion': 'SO:0000159',  # deletion
            'Inversion': 'SO:1000036',  # inversion
            'Not Applicable': 'SO:0001060',  # sequence variant - correct?
            'Not Specified': 'SO:0001060',  # sequence variant - correct?
            'Nucleotide repeat expansion': 'SO:0000667',  # insertion - correct?
            'Nucleotide substitutions': 'SO:1000002',  # substitution - Correct? Or another term indicating more than one?
            'Other': 'SO:0001060',  # sequence variant - correct?
            'Single point mutation': 'SO:1000008',  # point_mutation
            'Translocation': 'SO:0000199',  # translocation
            'Transposon insertion': 'SO:0000101',  # transposable_element
            'Undefined': 'SO:0001060',  # sequence variant - correct?
            'Viral insertion': 'SO:0000667',  # insertion - correct?
            'wild type': 'SO:0000817'  # wild type
        }
        if (sequence_alteration_type.strip() in type_map):
            type = type_map.get(sequence_alteration_type.strip())
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(allele_type)
            # print("Mapped: ", sequence_alteration_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Sequence Alteration Type (", sequence_alteration_type, ") not mapped")

        return type

    def _map_zygosity(self, zygosity):
        type = None
        type_map = {
            'Heterozygous': 'GENO:0000135',
            'Hemizygous Y-linked': 'GENO:0000604',
            'Heteroplasmic': 'GENO:0000603',
            'Homozygous': 'GENO:0000136',
            'Homoplasmic': 'GENO:0000602',
            'Hemizygous Insertion': 'GENO:0000606',
            'Hemizygous Deletion': 'GENO:0000606',  # hemizygous insertion
            #NOTE: GENO:0000606 is  'hemizygous insertion' but is used for the general 'hemizgous' in the Genotype.py file.
            'Hemizygous X-linked': 'GENO:0000605',
            'Indeterminate': 'GENO:0000137'
        }
        if (zygosity.strip() in type_map):
            type = type_map.get(zygosity)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Allele Type (", zygosity, ") not mapped")

        return type

    def _map_marker_type(self, marker_type):
        type = None
        type_map = {
            'Complex/Cluster/Region': 'SO:0000001',  # region. Something more specific available? #fixme
            'Transgene': 'SO:0000902',  # transgene
            'Gene': 'SO:0000704',  # gene
            'QTL': 'SO:0000771',  # QTL
            'DNA Segment': 'SO:0000110',  # sequence_feature. sequence_motif=SO:0001683? region=SO:0000001
            'Pseudogene': 'SO:0000336',  # pseudogene
            'Cytogenetic Marker': 'SO:0001645',  # genetic_marker?   #fixme
            'Other Genome Feature': 'SO:0000110',  # sequence_feature. Or sequence_motif=SO:0001683?
            'BAC/YAC end': 'SO:0000150',  # BAC_end: SO:0000999, YAC_end: SO:00011498; using parent term
        }
        if (marker_type.strip() in type_map):
            type = type_map.get(marker_type)
        else:
            print("ERROR: Marker Type (", marker_type, ") not mapped")

        return type

    def _map_taxon(self, taxon_name):
        type = None
        type_map = {
            'Bos taurus': 'NCBITaxon:9913',
            'Canis familiaris': 'NCBITaxon:9615',
            'Capra hircus': 'NCBITaxon:9925',
            'Cavia porcellus': 'NCBITaxon:10141',
            'Cricetulus griseus': 'NCBITaxon:10029',
            'Danio rerio': 'NCBITaxon:7955',
            'Equus caballus': 'NCBITaxon:9796',
            'Felis catus': 'NCBITaxon:9685',
            'Gallus gallus': 'NCBITaxon:9031',
            'Gorilla gorilla': 'NCBITaxon:9593',
            'Homo sapiens': 'NCBITaxon:9606',
            'Macaca mulatta': 'NCBITaxon:9544',
            'Macropus eugenii': 'NCBITaxon:9315',
            'Mesocricetus auratus': 'NCBITaxon:10036',
            'Microcebus murinus': 'NCBITaxon:30608',
            'Mus musculus/domesticus': 'NCBITaxon:10090',  # 10090=Mus musculus, 10092=Mus musculus domesticus
            'Ornithorhynchus anatinus': 'NCBITaxon:9258',
            'Oryctolagus cuniculus': 'NCBITaxon:9986',
            'Ovis aries': 'NCBITaxon:9940',
            'Pan troglodytes': 'NCBITaxon:9598',
            'Pongo pygmaeus': 'NCBITaxon:9600',
            'Rattus norvegicus': 'NCBITaxon:10116',
            'Sus scrofa domestica L.': 'NCBITaxon:9823',  # 9823=Sus scrofa, 9825=Sus scrofa domestica
            'Xenopus (Silurana) tropicalis': 'NCBITaxon:8364',
        }
        if (taxon_name.strip() in type_map):
            type = type_map.get(taxon_name)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Taxon Name (", taxon_name, ") not mapped")

        return type

    def _map_evidence_id(self, evidence_code):
        #TODO a default evidence code???  what should it be?
        type = None
        type_map = {
            'EXP': 'ECO:0000006',
            'IBA': 'ECO:0000318',
            'IC': 'ECO:0000001',
            'IDA': 'ECO:0000314',
            'IEA': 'ECO:0000501',
            'IEP': 'ECO:0000008',
            'IGI': 'ECO:0000316',
            'IKR': 'ECO:0000320',
            'IMP': 'ECO:0000315',
            'IPI': 'ECO:0000353',
            'ISA': 'ECO:0000200',
            'ISM': 'ECO:0000202',
            'ISO': 'ECO:0000201',
            'ISS': 'ECO:0000250',
            'NAS': 'ECO:0000303',
            'ND': 'ECO:0000035',
            'RCA': 'ECO:0000245',
            'TAS': 'ECO:0000304'
        }
        if (evidence_code.strip() in type_map):
            type = type_map.get(evidence_code)
        else:
            print("ERROR: Evidence code (", evidence_code, ") not mapped")

        return type

    def _map_evidence_label(self, evidence_code):
        type = None
        type_map = {
            'EXP': 'experimental evidence',
            'IBA': 'biological aspect of ancestor evidence used in manual assertion',
            'IC': 'inference from background scientific knowledge',
            'IDA': 'direct assay evidence used in manual assertion',
            'IEA': 'evidence used in automatic assertion',
            'IEP': 'expression pattern evidence',
            'IGI': 'genetic interaction evidence used in manual assertion',
            'IKR': 'phylogenetic determination of loss of key residues evidence used in manual assertion',
            'IMP': 'mutant phenotype evidence used in manual assertion',
            'IPI': 'physical interaction evidence used in manual assertion',
            'ISA': 'sequence alignment evidence',
            'ISM': 'match to sequence model evidence',
            'ISO': 'sequence orthology evidence',
            'ISS': 'sequence similarity evidence used in manual assertion',
            'NAS': 'non-traceable author statement used in manual assertion',
            'ND': 'no biological data found',
            'RCA': 'computational combinatorial evidence used in manual assertion',
            'TAS': 'traceable author statement used in manual assertion'
        }
        if (evidence_code.strip() in type_map):
            type = type_map.get(evidence_code)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Taxon Name (", evidence_code, ") not mapped")

        return type


    def _makeInternalIdentifier(self,prefix,key):
        '''
        This is a special MGI-to-MONARCH-ism.  MGI tables have unique keys that we use here, but don't want
        to necessarily re-distribute those internal identifiers.  Therefore, we make them into keys in a consistent
        way here.
        :param prefix: the object type to prefix the key with, since the numbers themselves are not unique across tables
        :param key: the number
        :return:
        '''

        return '_'+prefix+'key'+key

    def _querysparql(self):

        #load the graph
        vg = Graph()
        vg.parse(self.outfile, format="turtle")

        qres = g.query(
            """SELECT DISTINCT ?aname ?bname
                WHERE {
                    ?a foaf:knows ?b .
                    ?a foaf:name ?aname .
                    ?b foaf:name ?bname .
                }""")

        for row in qres:
            print("%s knows %s" % row)

        return