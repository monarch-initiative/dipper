from utils.CurieUtil import CurieUtil
from utils.GraphUtils import GraphUtils

__author__ = 'nicole'

from rdflib import Namespace
from rdflib.namespace import FOAF, DC, RDFS, OWL

import urllib, csv, os, time
from urllib import request
from datetime import datetime
from stat import *
import hashlib
import subprocess
from subprocess import check_call
from zipfile import ZipFile
import re

from models.InteractionAssoc import InteractionAssoc
from sources.Source import Source
from models.Dataset import Dataset

core_bindings = {'dc': DC, 'foaf': FOAF, 'rdfs': RDFS}



class BioGrid(Source):
    '''
    Biogrid interaction data
    '''

    files = {
        'interactions' : {'file': 'interactions.mitab.zip',
                           'url' : 'http://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.mitab.zip'
        },
        'identifiers' : {'file' : 'identifiers.tab.zip',
                          'url' : 'http://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip'
        }
    }

    curie_map = {
        'NCBIGene' : 'http://www.ncbi.nlm.nih.gov/gene/',
        'BIOGRID' : 'http://thebiogrid.org/',
        'WormBase' : 'http://identifiers.org/WormBase:',
        'MGI' : 'http://www.informatics.jax.org/accession/MGI:',
        'ENSEMBL' : 'http://identifiers.org/ENSEMBL:'
    }


    def __init__(self,args=[]):
        Source.__init__(self, 'biogrid')

        self.curie_map.update(InteractionAssoc.curie_map)

        self.load_bindings()

        self.dataset = Dataset('biogrid', 'The BioGrid', 'http://thebiogrid.org/')

        #data-source specific warnings (will be removed when issues are cleared)
        print("WARN: several MI experimental codes do not exactly map to ECO; using approximations.")
        return

    def fetch(self):
        '''
        :return: None
        '''

        for f in self.files.keys():
            file = self.files.get(f)
            self.fetch_from_url(file['url'],('/').join((self.rawdir,file['file'])))
            self.dataset.setFileAccessUrl(file['url'])

            st = os.stat(('/').join((self.rawdir,file['file'])))

        filedate=datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d")

        #the version number is encoded in the filename in the zip.
        #for example, the interactions file may unzip to BIOGRID-ALL-3.2.119.mitab.txt,
        #where the version number is 3.2.119
        f=('/').join((self.rawdir,self.files['interactions']['file']))
        with ZipFile(f, 'r') as myzip:
            flist=myzip.namelist()
            #assume that the first entry is the item
            fname=flist[0]
            #get the version from the filename
            version = re.match('BIOGRID-ALL-(\d+\.\d+\.\d+)\.mitab.txt',fname)
        myzip.close()

        self.dataset.setVersion(filedate,str(version.groups()[0]))

        return

    def parse(self,limit=None):
        '''
        abstract method to parse all data from an external resource, that was fetched in
        fetch()
        this should be overridden by subclasses
        :return: None
        '''

        #TODO make each of these items an option... we may want to process them separately
        #self._get_interactions(limit)

        self._get_identifiers(limit)

        self.load_bindings()


        ##### Write it out #####
        filewriter = open(self.outfile, 'w')
        self.load_bindings()
        print("Finished parsing files. Writing turtle to",self.outfile)
        print(self.graph.serialize(format="turtle").decode(),file=filewriter)
        filewriter.close()


        filewriter = open(self.datasetfile,'w')
        print(self.dataset.getGraph().serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        print("Wrote", len(self.graph), "nodes")


        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.curie_map.keys():
            v=self.curie_map[k]
            self.graph.bind(k, Namespace(v))
        return


    def verify(self):
        '''
        abstract method to verify the integrity of the data fetched and turned into triples
        this should be overridden by tests in subclasses
        :return: True if all tests pass
        '''
        status = False
        self._verify(self.outfile)
        status = self._verifyowl(self.outfile)

        return status

    def _get_interactions(self,limit):
        line_counter = 0
        f=('/').join((self.rawdir,self.files['interactions']['file']))
        myzip = ZipFile(f,'r')
        #assume that the first entry is the item
        fname=myzip.namelist()[0]
        matchcounter=0
        with myzip.open(fname,'r') as csvfile:
            for line in csvfile:
                #skip comment lines
                if (re.match('^#',line.decode())):
                    print("Skipping header line")
                    continue
                line_counter += 1
                line=line.decode().strip()
                #print(line)
                (interactor_a, interactor_b, alt_ids_a, alt_ids_b,aliases_a, aliases_b,
                 detection_method, pub_author,pub_id,
                 taxid_a, taxid_b, interaction_type,
                 source_db, interaction_id, confidence_val) = line.split('\t')

                #TODO remove these filters, or parameterize them
                #uncomment the following codeblock if you want to filter based on taxid
                #taxids = [9606,10090,10116,7227,7955,6239,8355]  #our favorite animals
                taxids = [9606] #human
                if (not (taxids.__contains__(int(re.sub('taxid:','', taxid_a.rstrip()))) or
                    taxids.__contains__(int(re.sub('taxid:','', taxid_b.rstrip()))) )):
                    continue
                else:
                    matchcounter += 1


                #TODO proper testing/catching of these search/match methods
                #get the actual gene ids, typically formated like: gene/locuslink:351|BIOGRID:106848
                gene_a='NCBIGene:'+re.search('locuslink\:(\d+)\|',interactor_a).groups()[0]
                gene_b='NCBIGene:'+re.search('locuslink\:(\d+)\|',interactor_b).groups()[0]

                #get the interaction type
                #psi-mi:"MI:0407"(direct interaction)
                int_type=re.search('MI:\d+',interaction_type).group()
                rel=self._map_MI_to_RO(int_type)

                #scrub pubmed-->PMID prefix
                pub_id = re.sub('pubmed','PMID',pub_id)
                #remove bogus whitespace
                pub_id = pub_id.strip()

                #get the method, and convert to evidence code
                det_code=re.search('MI:\d+',detection_method).group()
                evidence=self._map_MI_to_ECO(det_code)

                #note that the interaction_id is some kind of internal biogrid identifier that does not
                #map to a public URI.  we will construct a monarch identifier from this
                assoc_id = self.make_id(interaction_id)

                assoc = InteractionAssoc(assoc_id,gene_a,gene_b,pub_id,evidence,self.curie_map)
                assoc.setRelationship(rel)
                assoc.loadObjectProperties(self.graph)
                assoc.addInteractionAssociationToGraph(self.graph)
                if (limit is not None and line_counter > limit):
                    break

        myzip.close()
        #print("INFO: found",str(matchcounter),"matching rows")  #for testing

        return

    def _get_identifiers(self,limit):
        '''
        This will process the id mapping file provided by Biogrid.
        The file has a very large header, which we scan past, then pull the identifiers, and make
        equivalence axioms
        :param limit:
        :return:
        '''
        line_counter = 0
        f=('/').join((self.rawdir,self.files['identifiers']['file']))
        myzip = ZipFile(f,'r')
        #assume that the first entry is the item
        fname=myzip.namelist()[0]
        matchcounter=0
        foundheader=False
        with myzip.open(fname,'r') as csvfile:
            for line in csvfile:
                #skip header lines
                if (not foundheader):
                    if (re.match('BIOGRID_ID',line.decode())):
                        foundheader=True
                    continue

                line=line.decode().strip()
                #BIOGRID_ID	IDENTIFIER_VALUE	IDENTIFIER_TYPE	ORGANISM_OFFICIAL_NAME
                #1	814566	ENTREZ_GENE	Arabidopsis thaliana
                (biogrid_num,id_num,id_type,organism_label) = line.split('\t')

                #for each one of these, create the node and add equivalent classes
                cu = CurieUtil(self.curie_map)
                gu = GraphUtils(self.curie_map)
                g = self.graph
                biogrid_id='BIOGRID:'+biogrid_num
                prefix = self._map_idtype_to_prefix(id_type)


                #geneidtypefilters='NCBIGene,OMIM,MGI,FlyBase,ZFIN,MGI,HGNC,WormBase,XenBase,ENSEMBL,miRBase'.split(',')
                #TODO make these filters available as commandline options
                geneidtypefilters='NCBIGene,MGI,ENSEMBL'.split(',')
                #speciesfilters='Homo sapiens,Mus musculus,Drosophila melanogaster,Danio rerio,Caenorhabditis elegans,Xenopus laevis'.split(',')
                speciesfilters='Homo sapiens,Mus musculus'.split(',')
                #proteinidtypefilters='HPRD,Swiss-Prot,NCBIProtein'
                if ((speciesfilters is not None) and (organism_label.strip() in speciesfilters)):
                    line_counter += 1
                    if ((geneidtypefilters is not None) and (prefix in geneidtypefilters)):
                        mapped_id=(':').join((prefix,id_num))
                        gu.addEquivalentClass(g,biogrid_id,mapped_id)
                    elif (id_type == 'OFFICIAL_SYMBOL'):  #this symbol will only get attached to the biogrid class
                        gu.addClassToGraph(g,biogrid_id,id_num)
                    #elif (id_type == 'SYNONYM'):
                    #    gu.addSynonym(g,biogrid_id,id_num)  #FIXME - i am not sure these are synonyms, altids?

                if (limit is not None and line_counter > limit):
                    break

        return


    def _map_MI_to_RO(self,mi_id):
        mi_ro_map = {
            'MI:0403' : InteractionAssoc.relationships['colocalizes_with'], #colocalization
            'MI:0407' : InteractionAssoc.relationships['interacts_with'], #direct interaction
            'MI:0794' : InteractionAssoc.relationships['genetically_interacts_with'], #synthetic genetic interaction defined by inequality
            'MI:0796' : InteractionAssoc.relationships['genetically_interacts_with'], #suppressive genetic interaction defined by inequality
            'MI:0799' : InteractionAssoc.relationships['genetically_interacts_with'], #additive genetic interaction defined by inequality
            'MI:0914' : InteractionAssoc.relationships['interacts_with'], #association
            'MI:0915' : InteractionAssoc.relationships['interacts_with'] #physical association
        }

        ro_id = InteractionAssoc.relationships['interacts_with']  #default
        if (mi_id in mi_ro_map):
            ro_id = mi_ro_map.get(mi_id)

        return ro_id

    def _map_MI_to_ECO(self,mi_id):
        eco_id = 'ECO:0000006' #default to experimental evidence
        mi_to_eco_map = {
            'MI:0018' : 'ECO:0000068', #yeast two-hybrid
            'MI:0004' : 'ECO:0000079', #affinity chromatography
            'MI:0047' : 'ECO:0000076', #far western blotting
            'MI:0055' : 'ECO:0000021', #should be FRET, but using physical_interaction FIXME
            'MI:0090' : 'ECO:0000012', #desired: protein complementation, using: functional complementation
            'MI:0096' : 'ECO:0000085', #desired: pull down, using: immunoprecipitation
            'MI:0114' : 'ECO:0000324', #desired: x-ray crystallography, using: imaging assay
            'MI:0254' : 'ECO:0000011', #desired: genetic interference, using: genetic interaction evidence
            'MI:0401' : 'ECO:0000172', #desired: biochemical, using: biochemical trait evidence
            'MI:0415' : 'ECO:0000005', #desired: enzymatic study, using: enzyme assay evidence
            'MI:0428' : 'ECO:0000324', #imaging
            'MI:0686' : 'ECO:0000006', #desired: unspecified, using: experimental evidence
            'MI:1313' : 'ECO:0000006'  #None?
        }
        if (mi_id in mi_to_eco_map):
            eco_id = mi_to_eco_map.get(mi_id)
        else:
            print("WARN: unmapped code",mi_id, ". Defaulting to experimental_evidence")

        return eco_id

    def _map_idtype_to_prefix(self,idtype):
        prefix = idtype
        idtype_to_prefix_map = {
            'XENBASE' : 'XenBase',
            'TREMBL' : 'TrEMBL',
            'MGI' : 'MGI',
            'REFSEQ_DNA_ACCESSION': 'RefSeqNA',
            'MAIZEGDB' : 'MaizeGDB',
            'BEEBASE' : 'BeeBase',
            'ENSEMBL' : 'ENSEMBL',
            'TAIR' : 'TAIR',
            'GENBANK_DNA_GI' : 'NCBIgi',
            'CGNC' : 'CGNC',
            'RGD' : 'RGD',
            'GENBANK_GENOMIC_DNA_GI' : 'NCBIgi',
            'SWISSPROT' : 'Swiss-Prot',
            'MIM' : 'OMIM',
            'FLYBASE' : 'FlyBase',
            'VEGA' : 'VEGA',
            'ANIMALQTLDB' : 'AQTLDB',
            'ENTREZ_GENE_ETG' : 'ETG',
            'HPRD' : 'HPRD',
            'APHIDBASE' : 'APHIDBASE',
            'GENBANK_PROTEIN_ACCESSION' : 'NCBIProtein',
            'ENTREZ_GENE' : 'NCBIGene',
            'SGD' : 'SGD',
            'GENBANK_GENOMIC_DNA_ACCESSION' : 'NCBIGenome',
            'BGD' : 'BGD',
            'WORMBASE' : 'WormBase',
            'ZFIN' : 'ZFIN',
            'DICTYBASE' : 'dictyBase',
            'ECOGENE' : 'ECOGENE',
            'BIOGRID' : 'BIOGRID',
            'GENBANK_DNA_ACCESSION' : 'NCBILocus',
            'VECTORBASE' : 'VectorBase',
            'MIRBASE' : 'miRBase',
            'IMGT/GENE-DB' : 'IGMT',
            'HGNC' : 'HGNC',
            'SYSTEMATIC_NAME' : None,
            'OFFICIAL_SYMBOL' : None,
            'REFSEQ_GENOMIC_DNA_ACCESSION' : 'NCBILocus',
            'GENBANK_PROTEIN_GI' : 'NCBIgi',
            'REFSEQ_PROTEIN_ACCESSION' : 'RefSeqProt',
            'SYNONYM' : None,
            'GRID_LEGACY' : None
        }
        if (idtype in idtype_to_prefix_map):
            prefix = idtype_to_prefix_map.get(idtype)
        else:
            print("WARN: unmapped prefix",prefix)

        return prefix
