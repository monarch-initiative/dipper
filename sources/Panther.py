__author__ = 'nicole'

from rdflib.namespace import FOAF, DC, RDFS, OWL

import tarfile

import re

from models.Assoc import Assoc
from sources.Source import Source
from models.OrthologyAssoc import OrthologyAssoc
from models.Dataset import Dataset


class Panther(Source):
    '''
    The pairwise orthology calls from Panther DB: http://pantherdb.org/ encompass 22 species,
    from the RefGenome and HCOP projects.
    Here, we map the orthology classes to RO homology relationships
    This resource may be extended in the future with additional species.

    This currently makes a graph of orthologous relationships between genes, with the assumption that
    gene metadata (labels, equivalent ids) are provided from other sources.

    Gene families are nominally created from the orthology files, though these are incomplete with
    no hierarchical (subfamily) information.  This will get updated from the HMM files in the future.

    Note that there is a fair amount of identifier cleanup performed, to align with standard CURIE prefixes.
    '''

    files = {
        'refgenome' : {'file': 'RefGenomeOrthologs.tar.gz',
                           'url' : 'ftp://ftp.pantherdb.org/ortholog/current_release/RefGenomeOrthologs.tar.gz'
        },
        'hcop' : {'file' : 'orthologs_HCOP.tar.gz',
                          'url' : 'ftp://ftp.pantherdb.org/ortholog/current_release/orthologs_HCOP.tar.gz'
        }
    }

    namespaces = {
        'NCBIGene' : 'http://www.ncbi.nlm.nih.gov/gene/',
        'BIOGRID' : 'http://thebiogrid.org/',
        'ENSEMBL' : 'http://identifiers.org/ENSEMBL:',
        'MGI': 'http://www.informatics.jax.org/accession/MGI:',  #All MGI ids are genes in this case
        'UniProtKB' : 'http://identifiers.org/UniProt:',
        'RGD' : 'http://rgd.mcw.edu/rgdweb/report/gene/main.html?id=',
        'ZFIN' : 'http://zfin.org/',
        'dictyBase' : 'http://dictybase.org/gene/',
        'TAIR' : 'http://identifiers.org/TAIR:',
        'FlyBase' : 'http://identifiers.org/FB:',
        'PomBase' : 'http://identifiers.org/PomBase:',
        'WormBase' : 'http://identifiers.org/WormBase:',
        'SGD' : 'http://identifiers.org/SGD:',
        'PANTHER' : 'http://www.pantherdb.org/panther/family.do?clsAccession=',
        'DATA' : 'http://purl.obolibrary.org/DATA_',
    }


    def __init__(self,args=[]):
        Source.__init__(self, 'panther')

        self.namespaces.update(Assoc().get_namespaces())

        self.load_bindings()

        self.dataset = Dataset('panther', 'Protein ANalysis THrough Evolutionary Relationships', 'http://pantherdb.org/')

        #data-source specific warnings (will be removed when issues are cleared)


        return

    def fetch(self):
        '''
        :return: None
        '''

        self.get_files()
        #TODO the version number is tricky to get...we can't get it from redirects of the url
        #TODO use the remote timestamp of the file?

        return

    def parse(self,limit=None):
        '''
        abstract method to parse all data from an external resource, that was fetched in
        fetch()
        this should be overridden by subclasses
        :return: None
        '''

        self._get_orthologs(limit)

        self.load_bindings()


        ##### Write it out #####
        filewriter = open(self.outfile, 'w')
        self.load_bindings()
        print("INFO: Finished parsing files. Writing turtle to",self.outfile)
        print(self.graph.serialize(format="turtle").decode(),file=filewriter)
        filewriter.close()


        filewriter = open(self.datasetfile,'w')
        print(self.dataset.getGraph().serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        print("INFO: Wrote", len(self.graph), "nodes")


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

    def _get_orthologs(self,limit):
        '''
        will process each of the specified pairwise orthology files, creating orthology associations
        based on the specified orthology code.
        this currently assumes that each of the orthology files is identically formatted.

        there is also a nominal amount of identifier re-formatting:
        MGI:MGI --> MGI
        Ensembl --> ENSEMBL
        :param limit:
        :return:
        '''
        line_counter = 0
        for k in self.files.keys():
            f=('/').join((self.rawdir,self.files[k]['file']))
            matchcounter=0
            mytar = tarfile.open(f,'r:gz')

            #assume that the first entry is the item
            fname=mytar.getmembers()[0]
            print("INFO: Parsing",fname)
            with mytar.extractfile(fname) as csvfile:
                for line in csvfile:
                    #skip comment lines
                    if (re.match('^#',line.decode())):
                        print("INFO: Skipping header line")
                        continue
                    line_counter += 1
                    line=line.decode().strip()
                    #print(line)

                    # HUMAN|Ensembl=ENSG00000184730|UniProtKB=Q0VD83	MOUSE|MGI=MGI=2176230|UniProtKB=Q8VBT6	LDO	Euarchontoglires	PTHR15964

                    (a, b, orthology_class, ancestor_taxon,panther_id) = line.split('\t')
                    (species_a,gene_a,protein_a) = a.split('|')
                    (species_b,gene_b,protein_b) = b.split('|')

                    gene_a = re.sub('=',':',gene_a)
                    gene_b = re.sub('=',':',gene_b)

                    gene_a = self._clean_up_gene_id(gene_a,species_a)
                    gene_b = self._clean_up_gene_id(gene_b,species_b)

                    taxon_a = self._map_taxon_abbr_to_id(species_a)
                    taxon_b = self._map_taxon_abbr_to_id(species_b)

                    #TODO remove these filters, or parameterize them
                    ###uncomment the following code block if you want to filter based on taxid
                    taxids = [9606,10090,10116,7227,7955,6239,8355]  #our favorite animals
                    #taxids = [9606,10090] #human/mouse only
                    #taxids = [9606] #human only
                    #retain only those orthologous relationships to genes in the specified taxids
                    #using AND will get you only those associations where gene1 AND gene2 are in the taxid list (most-filter)
                    #using OR will get you any associations where gene1 OR gene2 are in the taxid list (some-filter)
                    if (not (taxids.__contains__(int(re.sub('NCBITaxon:','', taxon_a.rstrip()))) and
                        taxids.__contains__(int(re.sub('NCBITaxon:','', taxon_b.rstrip()))) )):
                        continue
                    else:
                        matchcounter += 1
                    ###end code block for filtering on taxon

                    rel=self._map_orthology_code_to_RO(orthology_class)

                    evidence = 'ECO:0000080'   #phylogenetic evidence

                    #note that the panther_id references a group of orthologs, and is not 1:1 with the rest
                    assoc_id = self.make_id(('').join((panther_id,species_a,gene_a,protein_a,species_b,gene_b,protein_b,orthology_class)))

                    assoc = OrthologyAssoc(assoc_id,gene_a,gene_b,None,evidence,self.namespaces)
                    assoc.setRelationship(rel)
                    assoc.loadObjectProperties(self.graph)
                    assoc.addAssociationToGraph(self.graph)

                    #note this is incomplete... it won't construct the full family hierarchy
                    assoc.addGeneFamilyToGraph(self.graph,(':').join(('PANTHER:',panther_id)))

                    if (limit is not None and line_counter > limit):
                        break

            print("INFO: finished processing",f)

        return


    def _map_taxon_abbr_to_id(self,ptax):
        '''
        Will map the panther-specific taxon abbreviations to NCBI taxon numbers
        :param ptax:
        :return: NCBITaxon id
        '''
        taxid = None
        ptax_to_taxid_map = {
            'HUMAN' : 9606,
            'SCHPO' : 4896,
            'ARATH' : 3702,
            'MOUSE' : 10090,
            'DANRE' : 7955,
            'YEAST' : 4932,
            'RAT' : 10116,
            'CAEEL' : 6239,
            'DICDI' : 44689,
            'CHICK' : 9031,
            'DROME' : 7227,
            'PIG' : 9823,
            'HORSE' : 9796,
            'CANFA' : 9615,
            'XENTR' : 8364,
            'TAKRU' : 31033,
            'MACMU' : 9544,
            'ORNAN' : 9258,
            'PANTR' : 9598,
            'ANOCA' : 28377,
            'MONDO' : 13616,
            'BOVIN' : 9913
        }
        if (ptax in ptax_to_taxid_map):
            taxid = (':').join(('NCBITaxon',str(ptax_to_taxid_map.get(ptax))))
        else:
            print("ERROR: unmapped taxon code",ptax)

        return taxid

    def _map_orthology_code_to_RO(self,ortho):
        '''
        Map the panther-specific orthology code (P,O,LDO,X,LDX) to relationship-ontology
        identifiers.
        :param ortho: orthology code
        :return: RO identifier
        '''
        ro_id = OrthologyAssoc.relationships['orthologous'] #in orthology relationship with
        ortho_to_ro_map = {
            'P' : OrthologyAssoc.relationships['paralogous'],
            'O' : OrthologyAssoc.relationships['orthologous'],
            'LDO': OrthologyAssoc.relationships['least_diverged_orthologous'],
            'X': OrthologyAssoc.relationships['xenologous'],
            'LDX': OrthologyAssoc.relationships['xenologous']
        }

        if (ortho in ortho_to_ro_map):
            ro_id = ortho_to_ro_map.get(ortho)
        else:
            print("WARN: unmapped orthology code",ortho,". Defaulting to 'orthology'.")

        return ro_id

    def _clean_up_gene_id(self,id,sp):
        '''
        A series of identifier rewriting to conform with standard gene identifiers.
        :param id:
        :param sp:
        :return:
        '''
        #special case for MGI
        id = re.sub('MGI:MGI:','MGI:',id)

        #rewrite Ensembl --> ENSEMBL
        id = re.sub('Ensembl','ENSEMBL',id)

        #rewrite Gene:CELE --> WormBase  these are old-school cosmid identifier
        id = re.sub('Gene:CELE','WormBase:',id)
        if (sp == 'CAEEL'):
            if (re.match('(Gene|ENSEMBLGenome):\w+\\.\d+',id)):
                id = re.sub('(?:Gene|ENSEMBLGenome):(\w+\\.\d+)','WormBase:\\1',id)

        #rewrite GeneID --> NCBIGene
        id = re.sub('GeneID','NCBIGene',id)

        #rewrite Gene:Dmel --> FlyBase
        id = re.sub('Gene:Dmel_','FlyBase:',id)
        #rewrite Gene:CG --> FlyBase:CG
        id = re.sub('Gene:CG','FlyBase:CG',id)

        #rewrite ENSEMBLGenome:FBgn --> FlyBase:FBgn
        id = re.sub('ENSEMBLGenome:FBgn','FlyBase:FBgn',id)

        #rewrite Gene:<ensembl ids> --> ENSEMBL:<id>
        id = re.sub('Gene:ENS','ENSEMBL:ENS',id)

        if (re.match('Gene:',id)):
            print("WARN: Found something I don't know how to fix ( species",sp,"):",id)

        return id