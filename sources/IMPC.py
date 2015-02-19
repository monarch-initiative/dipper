import csv
import gzip
import re
import logging

from rdflib import Literal
from rdflib.namespace import RDFS, RDF
from rdflib import URIRef

from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from models.G2PAssoc import G2PAssoc
from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils
from conf import curie_map


logger = logging.getLogger(__name__)


class IMPC(Source):

    files = {
        'impc': {'file': 'IMPC_genotype_phenotype.csv.gz',
                 'url': 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/IMPC_genotype_phenotype.csv.gz'},
        'euro': {'file': 'EuroPhenome_genotype_phenotype.csv.gz',
                 'url': 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/EuroPhenome_genotype_phenotype.csv.gz'},
        'mgd': {'file': 'MGP_genotype_phenotype.csv.gz',
                'url': 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/MGP_genotype_phenotype.csv.gz'},
        'checksum': {'file': 'checksum.md5',
                     'url': 'ftp://ftp.ebi.ac.uk/pub/databases/impc/latest/csv/checksum.md5'},
    }


    relationship = {
        'is_mutant_of': 'GENO:0000440',
        'derives_from': 'RO:0001000',
        'has_alternate_part': 'GENO:0000382',
        'has_reference_part': 'GENO:0000385',
        'in_taxon': 'RO:0000216',
        'has_zygosity': 'GENO:0000608',
        'is_sequence_variant_instance_of': 'GENO:0000408',
        'is_reference_instance_of': 'GENO:0000610',
        'hasExactSynonym': 'OIO:hasExactSynonym',
        'has_disposition': 'GENO:0000208',
        'has_phenotype': 'RO:0002200',
        'has_part': 'BFO:0000051',
        'has_variant_part': 'GENO:0000382'
    }

    terms = {
        'variant_locus': 'GENO:0000002',
        'reference_locus': 'GENO:0000036',
        'sequence_alteration': 'SO:0001059',
        'variant_single_locus_complement': 'GENO:0000030',
        'allele': 'GENO:0000008',
        'intrinsic_genotype': 'GENO:0000000',
        'effective_genotype': 'GENO:0000525',
        'phenotype': 'MONARCH:phenotype',  # Is this correct? What about GENO:0000348 - phenotype? MONARCH:phenotype
        'evidence': 'MONARCH:evidence',
        'genomic_background': 'GENO:0000010',
        'genomic_variation_complement': 'GENO:0000009',
        'zygosity': 'GENO:0000133'
    }

    def __init__(self):
        Source.__init__(self, 'impc')

        #update the dataset object with details about this resource
        #TODO put this into a conf file?
        self.dataset = Dataset('impc', 'IMPC', 'http://www.mousephenotype.org')

        #source-specific warnings.  will be cleared when resolved.
        #print("WARN: we are filtering G2P on the wild-type environment data for now")

        return


    def fetch(self, is_dl_forced):
        self.get_files(is_dl_forced)
        if self.compare_checksums():
            logger.debug('Files have same checksum as reference')
        else:
            raise Exception('Reference checksums do not match disk')
        return

    def scrub(self):
        '''
        Perform various data-scrubbing on the raw data files prior to parsing.
        For this resource, this currently includes:

        :return: None
        '''
        # scrub file of the oddities where there are "\" instead of empty strings
        #pysed.replace("\\\\", '', ('/').join((self.rawdir,self.files['geno']['file'])))

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        if (limit is not None):
            print("Only parsing first", limit, "rows of each file")
        print("Parsing files...")
        # IMPC data is delivered in three separate csv files. Need to iterate over each process (genotype, G2P, etc.)
        # for each file, unless instead we merge all three files into one, then process once.


        self._process_genotype_features(('/').join((self.rawdir,self.files['impc']['file'])), self.outfile, self.graph, limit)
        self._process_g2p(('/').join((self.rawdir, self.files['impc']['file'])), self.outfile, self.graph, limit)

        #TODO: Check processing of the other two IMPC data files
        self._process_genotype_features(('/').join((self.rawdir,self.files['euro']['file'])), self.outfile, self.graph, limit)
        self._process_g2p(('/').join((self.rawdir, self.files['euro']['file'])), self.outfile, self.graph, limit)
        self._process_genotype_features(('/').join((self.rawdir,self.files['mgd']['file'])), self.outfile, self.graph, limit)
        self._process_g2p(('/').join((self.rawdir, self.files['mgd']['file'])), self.outfile, self.graph, limit)


        print("Finished parsing.")

        self.load_bindings()
        Assoc().loadObjectProperties(self.graph)

        print("Found", len(self.graph), "nodes")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):


        print("Processing Genotypes")
        gu = GraphUtils(curie_map.get())
        cu = CurieUtil(curie_map.get())

        line_counter = 0
        with gzip.open(raw, 'rt') as csvfile:
        #with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                line_counter += 1
                #print(row)

                (resource_name,phenotyping_center,colony_id,strain_name,strain_accession_id,marker_symbol,
                 marker_accession_id,allele_symbol,allele_accession_id,zygosity,sex,pipeline_name,pipeline_stable_id,
                 procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,mp_term_name,mp_term_id,p_value,
                 effect_size) = row

                # For IMPC, we put together genotype IDs from other parts (marker, allele, and strain IDs, with
                # adjustments based on zygosity).



                # Making genotype IDs from the various parts, can change later if desired.
                # Should these just be the genotype IDs, or should we use the self._makeInternalIdentifier()?
                # Or perhaps this, adjusted base on zygosity:
                # genotype_id = self.make_id((marker_accession_id+allele_accession_id+strain_accession_id))
                # OR, even more simplified:
                # genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))
                # If we use this even more simplified method for the genotype ID,
                # may want to move the zygosity warning to the zygosity processing portion.

                #TODO: better to swap out the zygosity for the zygosity_id in all make_id statements?
                zygosity_id = self._map_zygosity(zygosity)
                zygosity_name = self._map_zygosity_name(zygosity)

                genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))

                # Make the variant locus name/label
                if re.match(".*<.*>.*", allele_symbol):
                    variant_locus_name = allele_symbol
                else:
                    variant_locus_name = allele_symbol+'<'+allele_symbol+'>'

                variant_locus_id = self.make_id((marker_accession_id+allele_accession_id))
                variant_locus_type = self.terms['variant_locus']

                # Making VSLC labels from the various parts, can change later if desired.
                if zygosity == 'heterozygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*', '<+>', variant_locus_name)
                    #print(vslc_name)
                elif zygosity == 'homozygote':
                    vslc_name = variant_locus_name+'/'+variant_locus_name
                    #print(vslc_name)
                elif zygosity == 'hemizygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*', '<0>', variant_locus_name)
                    #print(vslc_name)
                elif zygosity == 'not_applicable':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*', '<?>', variant_locus_name)
                    #print(vslc_name)
                # Do we need to handle an unknown zygosity label in another fashion?
                # How do we want to handle the "not_applicable" zygosity labels?
                # Is the question mark allele the correct way?
                else:
                    print("INFO: found unknown zygosity :",zygosity)
                    break

                # Create the GVC labels
                # For IMPC, the GVC == VSLC
                gvc_name = vslc_name


                # Making genotype labels from the various parts, can change later if desired.
                genotype_name = vslc_name+'['+strain_name+']'


                # Create the genotype
                geno = Genotype(genotype_id, genotype_name)
                #print('ID='+genotype_id+'     LABEL='+genotype_name)



                #This is for handling any of the alleles that do not have an MGI ID or IMPC has not yet
                # updated their data for the allele with the MGI ID. The IDs look like NULL-<10-digit string>.

                if re.match("MGI:.*",allele_accession_id):
                    sequence_alteration_id = allele_accession_id
                else:
                    sequence_alteration_id = 'IMPC:'+allele_accession_id

                if re.match(".*<.*",allele_symbol):
                    sequence_alteration_name = re.sub('.*<','<',allele_symbol)
                else:
                    sequence_alteration_name = '<'+allele_symbol+'>'


                #print()
                #Add allele to genotype
                #FIXME: Is it correct to add the type to the allele, or should this be added to the sequence alteration?
                #Hijacking this to be the variant_locus instead of the allele type.
                geno.addAllele(variant_locus_id, variant_locus_name, variant_locus_type, None)

                #Hard coding gene_type as gene.
                gene_type = 'SO:0000704'# gene

                if re.match('MGI:.*',marker_accession_id):
                    gene_id = marker_accession_id
                else:
                    gene_id = 'IMPC:'+marker_accession_id


                #Add gene to genotype
                #need gene_id, gene_label, gene_type. No gene description in data.
                #All marker IDs have the MGI: prefix, but should we include an error message should a
                # future data set not be in that format?
                geno.addGene(gene_id, marker_symbol,gene_type)

                #FIXME: What about heterozygotes? Do we need to add wild type alleles? No identifier included.
                # Add allele to gene
                geno.addAlleleOfGene(variant_locus_id, gene_id)

                # genotype has_part allele
                #FIXME: What about heterozygotes?  Do we need to add wild type alleles? No identifier included.
                geno.addAlleleToGenotype(genotype_id, variant_locus_id)

                # need to make some attributes of this relationship for allele zygosity within the genotype
                #or do we just make the allele_complement here, based on the zygosity?
                #allele_in_gene_id=self.make_id(genotype_id+allele_id+zygosity)
                #allele has_disposition zygosity?
                #                g.add(())

                #TODO
                #The following parts are test code for creating the more complex parts of the genotype partonomy.
                #Will work on abstracting these code snippets to generalized classes for use by any resource.

                # Add the VSLC
                # Create the VSLC ID
                vslc_id = self.make_id((marker_accession_id+allele_accession_id+zygosity))
                # vslc is of type vslc
                self.graph.add((URIRef(cu.get_uri(vslc_id)),RDF['type'],URIRef(cu.get_uri(self.terms['variant_single_locus_complement']))))
                #vslc has label vslc_name
                self.graph.add((URIRef(cu.get_uri(vslc_id)),RDFS['label'],Literal(vslc_name)))
                #vslc has_alternate_part allele/variant_locus
                #FIXME: Matt's diagram has this one starred, perhaps still deciding on the has_part vs has_variant_part designation?
                self.graph.add((URIRef(cu.get_uri(vslc_id)),URIRef(cu.get_uri(self.relationship['has_variant_part'])),URIRef(cu.get_uri(variant_locus_id))))

                # Add the zygosity
                # Add zygosity as a class
                gu.addClassToGraph(self.graph,zygosity_id,zygosity_name)
                # Zygosity is of type zygosity #FIXME: Is this correct?
                self.graph.add((URIRef(cu.get_uri(zygosity_id)),RDF['type'],URIRef(cu.get_uri(self.terms['zygosity']))))
                # Add the zygosity to the VSLC
                # Need a connection based on GENO:0000608 has_zygosity
                self.graph.add(((URIRef(cu.get_uri(vslc_id)),URIRef(cu.get_uri(self.relationship['has_zygosity'])),URIRef(cu.get_uri(zygosity_id)))))


                # Link the VSLC to the allele based on zygosity
                # Need a connection based on GENO:0000382 has_alternate_part or
                # GENO:0000385 has_reference_part.
                #FIXME:What about for hemizygote/unknown?


                # Add the GVC
                # Create the GVC ID
                gvc_id = self.make_id(('GVC:'+marker_accession_id+allele_accession_id+zygosity_id))
                # gvc is of type gvc
                self.graph.add((URIRef(cu.get_uri(gvc_id)),RDF['type'],URIRef(cu.get_uri(self.terms['genomic_variation_complement']))))
                #gvc has label gvc
                self.graph.add((URIRef(cu.get_uri(gvc_id)),RDFS['label'],Literal(gvc_name)))
                #gvc has_variant_part vslc
                self.graph.add((URIRef(cu.get_uri(gvc_id)),URIRef(cu.get_uri(self.relationship['has_variant_part'])),URIRef(cu.get_uri(vslc_id))))
                # Add the GVC to the intrinsic genotype
                self.graph.add((URIRef(cu.get_uri(genotype_id)),URIRef(cu.get_uri(self.relationship['has_variant_part'])),URIRef(cu.get_uri(gvc_id))))


                # Add the effective genotype
                effective_genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id+sex))
                effective_genotype_label = genotype_name+'('+sex+')'
                 # effective_genotype is of type effective_genotype
                self.graph.add((URIRef(cu.get_uri(effective_genotype_id)),RDF['type'],URIRef(cu.get_uri(self.terms['effective_genotype']))))
                #effective_genotype has label effective_genotype_label
                self.graph.add((URIRef(cu.get_uri(effective_genotype_id)),RDFS['label'],Literal(effective_genotype_label)))
                # Add the intrinsic_genotype to the effective_genotype
                self.graph.add((URIRef(cu.get_uri(effective_genotype_id)),URIRef(cu.get_uri(self.relationship['has_alternate_part'])),URIRef(cu.get_uri(genotype_id))))

                # Add the genomic background
                 # create the genomic background id and name
                if re.match("MGI:.*",strain_accession_id):
                    genomic_background_id = strain_accession_id
                else:
                    genomic_background_id = 'IMPC:'+strain_accession_id
                    #FIXME: Will this resolve, or do we need a separate IMPCStrain:?
                genomic_background_name = strain_name
                # genomic_background_id is of type genomic_background
                self.graph.add((URIRef(cu.get_uri(genomic_background_id)),RDF['type'],URIRef(cu.get_uri(self.terms['genomic_background']))))
                # genomic_background_id has label genomic_background_name
                self.graph.add((URIRef(cu.get_uri(genomic_background_id)),RDFS['label'],Literal(genomic_background_name)))
                # intrinsic_genotype_id has_reference_part genomic_background_id
                self.graph.add((URIRef(cu.get_uri(genomic_background_id)),URIRef(cu.get_uri(self.relationship['has_reference_part'])),URIRef(cu.get_uri(genomic_background_id))))


                # Add the taxon as a class
                #FIXME: My cmap has it indicated as a class. Is that correct? EDIT: Think so, MB's map has it as a class
                # Will this be redundant with ingest of NCBITaxonomy?
                taxon_id = 'NCBITaxon:10090'
                taxon_name = 'Mus musculus'
                gu.addClassToGraph(self.graph,taxon_id,taxon_name)
                # Add the taxon to the genomic_background_id
                self.graph.add((URIRef(cu.get_uri(genomic_background_id)),URIRef(cu.get_uri(self.relationship['in_taxon'])),URIRef(cu.get_uri(taxon_id))))

                #Add the sequence_alteration and sequence_alteration_type
                #FIXME
                #IMPC contains targeted mutations with either gene traps, knockouts, insertion/intragenic deletions.
                #Currently hard-coding to the insertion type. Does this need to be adjusted?
                sequence_alteration_type_id = 'SO:0000667'  # insertion
                sequence_alteration_type_name = 'insertion'  # insertion
                 # sequence_alteration is of type sequence_alteration
                self.graph.add((URIRef(cu.get_uri(sequence_alteration_id)),RDF['type'],URIRef(cu.get_uri(self.terms['sequence_alteration']))))
                #sequence_alteration has label sequence_alteration_label
                self.graph.add((URIRef(cu.get_uri(sequence_alteration_id)),RDFS['label'],Literal(sequence_alteration_name)))
                # Add sequence_alteration_type as a class
                gu.addClassToGraph(self.graph,sequence_alteration_type_id,sequence_alteration_type_name)
                #sequence_alteration has type sequence_alteration_type
                self.graph.add((URIRef(cu.get_uri(sequence_alteration_id)),RDF['type'],URIRef(cu.get_uri(sequence_alteration_type_id))))
                # Add the sequence_alteration to the variant_locus
                self.graph.add((URIRef(cu.get_uri(variant_locus_id)),URIRef(cu.get_uri(self.relationship['has_variant_part'])),URIRef(cu.get_uri(sequence_alteration_id))))

                #TODO: sequence_alteration_type



                #add the specific genotype subgraph to the overall graph
                if (limit is not None and line_counter > limit):
                    break
                self.graph = geno.getGraph().__iadd__(self.graph)

        print("INFO: Done with genotypes")
        return



    def _process_g2p(self, raw, out, g, limit=None):
        '''
        This module currently filters for only wild-type environments, which clearly excludes application
        of morpholinos.  Very stringent filter.  To be updated at a later time.
        :param raw:
        :param out:
        :param g:
        :param limit:
        :return:
        '''

        #NOTE: No evidence provided in the data file, no environment data

        #QUESTION: A matter of efficiency, but is it better to process the genotype ang G2P sections separately,
        # or in a case like IMPC with single files, we could potentially process the G2P portions within the Genotype call.
        # Current implementation results in processing the files multiple times, which isn't terrible given that
        # the files for IMPC are not large files.


        print("Processing G2P")


        line_counter = 0
        with gzip.open(raw, 'rt') as csvfile:
        #with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            next(filereader, None)  # skip the header row
            for row in filereader:
                line_counter += 1
                #print(row)

                (resource_name,phenotyping_center,colony_id,strain_name,strain_accession_id,marker_symbol,
                 marker_accession_id,allele_symbol,allele_accession_id,zygosity,sex,pipeline_name,pipeline_stable_id,
                 procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,mp_term_name,mp_term_id,p_value,
                 effect_size) = row


                phenotype_id = mp_term_id

                #Not currently using the phenotype name, but add it somewhere as the label?
                phenotype_name = mp_term_name

                # Making genotype IDs from the various parts, can change later if desired.
                # Should these just be the genotype IDs, or should we use the self._makeInternalIdentifier()?
                genotype_id = self.make_id((marker_accession_id+allele_accession_id+zygosity+strain_accession_id))

                # Make the variant locus name/label
                if re.match(".*<.*>.*", allele_symbol):
                    variant_locus_name = allele_symbol
                else:
                    variant_locus_name = allele_symbol+'<'+allele_symbol+'>'

                # Making VSLC labels from the various parts, can change later if desired.
                if zygosity == 'heterozygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<+>',variant_locus_name)
                elif zygosity == 'homozygote':
                    vslc_name = variant_locus_name+'/'+variant_locus_name
                elif zygosity == 'hemizygote':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<0>',variant_locus_name)
                elif zygosity == 'not_applicable':
                    vslc_name = variant_locus_name+'/'+re.sub('<.*','<?>',variant_locus_name)
                else:
                    print("INFO: found unknown zygosity :",zygosity)
                    break

                # Create the GVC labels
                # For IMPC, the GVC == VSLC
                gvc_name = vslc_name

                # Making genotype labels from the various parts, can change later if desired.
                genotype_name = gvc_name+'['+strain_name+']'

                geno = Genotype(genotype_id, genotype_name)
                self.graph.__iadd__(geno.getGraph())

                #FIXME
                #In the NIF/DISCO view, we don't really have a publication id. We have a publication URL.
                #pub_url example: http://www.mousephenotype.org/data/genes/MGI:108077
                #There are duplicates, and will be duplicates of the genotype/phenotype/pub_id combination.
                #but could add the testing parameter?
                pub_id = 'IMPC:'+marker_accession_id

                #FIXME
                #No evidence code provided. NIF/DISCO view was hard coded to null. However,
                # isn't all of the IMPC data based on experimental evidence? Or is that too general?
                # Could use 'EXP': 'ECO:0000006', although the code below used in ZFIN might be more appropriate.
                eco_id = "ECO:0000059"  # experimental_phenotypic_evidence This was used in ZFIN

                #pub_id could be removed here as well.
                assoc_id = self.make_id((genotype_id+phenotype_id+pub_id))

                #removed pub_id.
                #assoc = G2PAssoc(assoc_id, genotype_id, phenotype_id, pub_id, eco_id)
                assoc = G2PAssoc(assoc_id, genotype_id, phenotype_id, None, eco_id)
                self.graph = assoc.addAssociationNodeToGraph(self.graph)

                #TODO: Add the additional phenotype data parts.
                #Will have to figure out exactly how to map these relationships, abstract to reusable functions.

                #In the G2PAssoc, the association is mapped between the genotype and the phenotype, but the
                # phenotype is not really created as a separate instance that we can then hang additional
                # phenotype parts on, correct? HOWEVER, while the phenotype term ID is not unique,
                # the phenotype as measured by a specific testing parameter in a specific effective genotype
                # would be unique. Are there multiple testing parameters for that combination?

                #Add phenotype_description_free_text
                #Format: 'Phenotype observed by '|| phenotyping_center||' in an '||procedure_name||' assay where '||testing_parameter_name||' was measured with an effect_size of '||effect_size||'. (p='||p_value||')' as phenotype_description_free_text,
                phenotype_description_free_text = 'Phenotype observed by '+phenotyping_center+' in an '+procedure_name+' assay where '+parameter_name+' was measured with an effect_size of '+effect_size+'. (p='+p_value+')'


                #Add phenotype description to phenotype
                #g.add((n, DC['description'], Literal(phenotype_description_free_text)))


                # Add the phenotype
                # Create the phenotype ID. Can this be repeated given that the phenotype terms are very common,
                # or do we need to create a unique phenotype id for each row?
                #vslc_id = self.make_id((marker_accession_id+allele_accession_id+zygosity))
                # vslc is of type vslc
                #self.graph.add((URIRef(cu.get_uri(vslc_id)),RDF['type'],URIRef(cu.get_uri(self.terms['variant_single_locus_complement']))))
                #vslc has label vslc_name
                #self.graph.add((URIRef(cu.get_uri(vslc_id)),RDFS['label'],Literal(vslc_name)))
                #vslc has_alternate_part allele/variant_locus

                #self.graph.add((URIRef(cu.get_uri(vslc_id)),URIRef(cu.get_uri(self.relationship['has_variant_part'])),URIRef(cu.get_uri(variant_locus_id))))




                if (limit is not None and line_counter > limit):
                    break



        return


    def _map_zygosity(self, zygosity):
        type = None
        type_map = {
            'heterozygote': 'GENO:0000135',
            'homozygote': 'GENO:0000136',
            'hemizygote': 'GENO:0000134',
            'not_applicable': 'GENO:0000137'
        }
        if (zygosity.strip() in type_map):
            type = type_map.get(zygosity)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Zygosity Type (", zygosity, ") not mapped")

        return type

    def _map_zygosity_name(self, zygosity):
        type = None
        type_map = {
            'heterozygote': 'heterozygous',
            'homozygote': 'homozygous',
            'hemizygote': 'hemizygous',
            'not_applicable': 'unspecified zygosity'
        }
        if (zygosity.strip() in type_map):
            type = type_map.get(zygosity)
            # type = 'http://purl.obolibrary.org/obo/' + type_map.get(zygosity)
        # print("Mapped: ", allele_type, "to", type)
        else:
            # TODO add logging
            print("ERROR: Zygosity Type (", zygosity, ") not mapped")

        return type

    def parse_checksum_file(self,file):
        """
        :param file
        :return dict
        """
        checksums = dict()
        file_path = '/'.join((self.rawdir, file))
        with open(file_path, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter=' ')
            for row in reader:
                (checksum, whitespace, file_name) = row
                checksums[checksum] = file_name

        return checksums

    def compare_checksums(self):
        """
        test to see if fetched file matches checksum from ebi
        :return: True or False
        """
        is_match = True
        reference_checksums = self.parse_checksum_file(
            self.files['checksum']['file'])
        for md5, file in reference_checksums.items():
            if self.get_file_md5(self.rawdir, file) != md5:
                is_match = False
                logger.warn('%s was not downloaded completely', file)
                return is_match

        return is_match

    def verify(self):
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        # verify some kind of relationship that should be in the file
        return status


