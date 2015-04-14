__author__ = 'nicole'

import tarfile
import re

from dipper.sources.Source import Source
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper import config, curie_map


class Panther(Source):
    """
    The pairwise orthology calls from Panther DB: http://pantherdb.org/ encompass 22 species,
    from the RefGenome and HCOP projects.
    Here, we map the orthology classes to RO homology relationships
    This resource may be extended in the future with additional species.

    This currently makes a graph of orthologous relationships between genes, with the assumption that
    gene metadata (labels, equivalent ids) are provided from other sources.

    Gene families are nominally created from the orthology files, though these are incomplete with
    no hierarchical (subfamily) information.  This will get updated from the HMM files in the future.

    Note that there is a fair amount of identifier cleanup performed, to align with our standard CURIE prefixes.

    The test graph of data is output based on configured "protein" identifiers in conf.json.

    """

    files = {
        'refgenome' : {'file': 'RefGenomeOrthologs.tar.gz',
                           'url' : 'ftp://ftp.pantherdb.org/ortholog/current_release/RefGenomeOrthologs.tar.gz'
        },
        'hcop' : {'file' : 'orthologs_HCOP.tar.gz',
                          'url' : 'ftp://ftp.pantherdb.org/ortholog/current_release/orthologs_HCOP.tar.gz'
        }
    }

    def __init__(self, tax_ids=None):
        super().__init__('panther')
        self.tax_ids = tax_ids
        self.load_bindings()

        self.dataset = Dataset('panther', 'Protein ANalysis THrough Evolutionary Relationships', 'http://pantherdb.org/')

        # Defaults
        if self.tax_ids is None:
            self.tax_ids = [9606, 10090]

        if (not (('test_ids' in config.get_config()) and ('protein' in config.get_config()['test_ids']))):
            print("WARN: not configured with gene test ids.")
        else:
            self.test_ids = config.get_config()['test_ids']['protein']



        #data-source specific warnings (will be removed when issues are cleared)

        return

    def fetch(self, is_dl_forced):
        """
        :return: None
        """

        self.get_files(is_dl_forced)
        #TODO the version number is tricky to get...we can't get it from redirects of the url
        #TODO use the remote timestamp of the file?

        return

    def parse(self,limit=None):
        """
        abstract method to parse all data from an external resource, that was fetched in
        fetch()
        :return: None
        """

        loops = [True]
        if not self.testOnly:
            loops = [True,False]

        for testMode in loops:
            self._get_orthologs(limit,testMode)

        self.load_bindings()

        print("INFO: Found", len(self.graph), "nodes")


        return


    def _get_orthologs(self,limit,testMode):
        """
        This will process each of the specified pairwise orthology files, creating orthology associations
        based on the specified orthology code.
        this currently assumes that each of the orthology files is identically formatted.
        relationships are made between genes here.

        there is also a nominal amount of identifier re-formatting:
        MGI:MGI --> MGI
        Ensembl --> ENSEMBL

        we skip any genes where we don't know how to map the gene identifiers.  for example,
        Gene:Huwe1 for RAT is not an identifier, so we skip any mappings to this identifier.  Often, the
        there are two entries for the same gene (base on equivalent Uniprot id), and so we are not
        actually losing any information.

        We presently have a hard-coded filter to select only orthology relationships where one of the pair
        is in our species of interest (Mouse and Human, for the moment).  This will be added as a
        configurable parameter in the future.

        Genes are also added to a grouping class defined with a PANTHER id.

        Triples:
        <gene1_id> RO:othologous <gene2_id>
        <assoc_id> :hasSubject <gene1_id>
        <assoc_id> :hasObject <gene2_id>
        <assoc_id> :hasPredicate <RO:orthologous>
        <assoc_id> dc:evidence ECO:phylogenetic_evidence

        <panther_id> a DATA:gene_family
        <panther_id> RO:has_member <gene1_id>
        <panther_id> RO:has_member <gene2_id>

        :param limit:
        :return:
        """
        print("INFO: getting orthologs")
        line_counter = 0

        if testMode:
            g = self.testgraph

        else:
            g = self.graph

        gu = GraphUtils(curie_map.get())

        for k in self.files.keys():
            f=('/').join((self.rawdir,self.files[k]['file']))
            matchcounter=0
            mytar = tarfile.open(f,'r:gz')

            #assume that the first entry is the item
            fname=mytar.getmembers()[0]
            print("INFO: Parsing",fname.name)
            with mytar.extractfile(fname) as csvfile:
                for line in csvfile:
                    #skip comment lines
                    if (re.match('^#',line.decode())):
                        print("INFO: Skipping header line")
                        continue
                    line_counter += 1
                    line=line.decode().strip()
                    #print(line)

                    #parse each row
                    # HUMAN|Ensembl=ENSG00000184730|UniProtKB=Q0VD83	MOUSE|MGI=MGI=2176230|UniProtKB=Q8VBT6	LDO	Euarchontoglires	PTHR15964
                    (a, b, orthology_class, ancestor_taxon,panther_id) = line.split('\t')
                    (species_a,gene_a,protein_a) = a.split('|')
                    (species_b,gene_b,protein_b) = b.split('|')

                    #skip the entries that don't have homolog relationships with the test ids
                    if testMode and not (re.sub('UniProtKB=','',protein_a) in self.test_ids or
                        re.sub('UniProtKB=','',protein_b) in self.test_ids):
                        continue


                    #map the taxon abbreviations to ncbi taxon ids
                    taxon_a = self._map_taxon_abbr_to_id(species_a)
                    taxon_b = self._map_taxon_abbr_to_id(species_b)

                    #TODO remove these filters, or parameterize them
                    ###uncomment the following code block if you want to filter based on taxid
                    #taxids = [9606,10090,10116,7227,7955,6239,8355]  #our favorite animals
                    #taxids = [9606] #human only
                    #retain only those orthologous relationships to genes in the specified taxids
                    #using AND will get you only those associations where gene1 AND gene2 are in the taxid list (most-filter)
                    #using OR will get you any associations where gene1 OR gene2 are in the taxid list (some-filter)
                    #print("INFO: restricting taxa to ids:",taxids)
                    if ((int(re.sub('NCBITaxon:','', taxon_a.rstrip())) not in self.tax_ids ) and
                        (int(re.sub('NCBITaxon:','', taxon_b.rstrip())) not in self.tax_ids ) ):
                        continue
                    else:
                        matchcounter += 1
                        if (limit is not None and matchcounter > limit):
                            break

                    ###end code block for filtering on taxon

                    #fix the gene identifiers
                    gene_a = re.sub('=',':',gene_a)
                    gene_b = re.sub('=',':',gene_b)

                    gene_a = self._clean_up_gene_id(gene_a,species_a)
                    gene_b = self._clean_up_gene_id(gene_b,species_b)

                    #a special case here; mostly some rat genes they use symbols instead of identifiers.  will skip
                    if (re.match('Gene:',gene_a) or re.match('Gene:',gene_b)):
                        continue

                    rel=self._map_orthology_code_to_RO(orthology_class)

                    evidence = 'ECO:0000080'   #phylogenetic evidence

                    #note that the panther_id references a group of orthologs, and is not 1:1 with the rest
                    assoc_id = self.make_id(('').join((panther_id,species_a,gene_a,protein_a,species_b,gene_b,protein_b,orthology_class)))

                    #add the association and relevant nodes to graph
                    assoc = OrthologyAssoc(assoc_id,gene_a,gene_b,None,evidence)
                    assoc.setRelationship(rel)
                    assoc.loadAllProperties(g)    #FIXME inefficient

                    #add genes to graph; assume labels will be taken care of elsewhere
                    gu.addClassToGraph(g,gene_a,None)
                    gu.addClassToGraph(g,gene_b,None)

                    assoc.addAssociationToGraph(g)

                    #note this is incomplete... it won't construct the full family hierarchy, just the top-grouping
                    assoc.addGeneFamilyToGraph(self.graph,(':').join(('PANTHER',panther_id)))

                    if (not testMode) and (limit is not None and line_counter > limit):
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
        o = OrthologyAssoc(None,None,None,None,None)  #this is sneaky, needs refactor

        ro_id = o.properties['orthologous'] #in orthology relationship with
        ortho_to_ro_map = {
            'P' : o.properties['paralogous'],
            'O' : o.properties['orthologous'],
            'LDO': o.properties['least_diverged_orthologous'],
            'X': o.properties['xenologous'],
            'LDX': o.properties['xenologous']
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