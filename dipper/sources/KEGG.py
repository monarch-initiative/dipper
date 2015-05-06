import csv
import logging
import re

from dipper.sources.Source import Source
from dipper.models.Dataset import Dataset
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.Genotype import Genotype
from dipper.models.OrthologyAssoc import OrthologyAssoc
from dipper.utils.GraphUtils import GraphUtils
from dipper import curie_map

logger = logging.getLogger(__name__)


class KEGG(Source):

    files = {
        'disease': {'file': 'disease',
                 'url': 'http://rest.genome.jp/list/disease'},
        'pathway': {'file': 'pathway',
                 'url': 'http://rest.genome.jp/list/pathway'},
        'hsa_genes': {'file': 'hsa_genes',
                 'url': 'http://rest.genome.jp/list/hsa'},
        'ortholog_classes': {'file': 'ortholog_classes',
                 'url': 'http://rest.genome.jp/list/orthology'},
        'disease_gene': {'file': 'disease_gene',
                 'url': 'http://rest.kegg.jp/link/disease/hsa'},
        'omim2disease': {'file': 'omim2disease',
                 'url': 'http://rest.genome.jp/link/disease/omim'},
        'omim2gene': {'file': 'omim2gene',
                 'url': 'http://rest.genome.jp/link/omim/hsa'},
        'ncbi': {'file': 'ncbi',
                 'url': 'http://rest.kegg.jp/conv/ncbi-geneid/hsa'},
        'hsa_gene2pathway': {'file': 'human_gene2pathway',
                 'url': 'http://rest.kegg.jp/link/pathway/hsa'},
        'hsa_orthologs': {'file': 'hsa_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/hsa'},
        'mmu_orthologs': {'file': 'mmu_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/mmu'},
        'rno_orthologs': {'file': 'rno_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/rno'},
        'dme_orthologs': {'file': 'dme_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/dme'},
        'dre_orthologs': {'file': 'dre_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/dre'},
        'cel_orthologs': {'file': 'cel_orthologs',
                 'url': 'http://rest.kegg.jp/link/orthology/cel'}
    }

    # I do not love putting these here; but I don't know where else to put them
    test_ids = {
        "pathway": ["path:map00010", "path:map00195", "path:map00100", "path:map00340"],
        "disease": ["ds:H00015", "ds:H00026", "ds:H00712", "ds:H00736"],
        "genes": ["hsa:100506275", "hsa:285958", "hsa:286410", "hsa:6387"],
        "orthology_classes": ["ko:K00010", "ko:K00027", "ko:K00042", "ko:K00088"]
    }

    def __init__(self):
        Source.__init__(self, 'kegg')

        # update the dataset object with details about this resource
        # TODO put this into a conf file?
        self.dataset = Dataset('kegg', 'KEGG', 'http://www.genome.jp/kegg/', None, None)

        # source-specific warnings.  will be cleared when resolved.

        return

    def fetch(self, is_dl_forced):
        self.get_files(is_dl_forced)
        #if self.compare_checksums():
            #logger.debug('Files have same checksum as reference')
        #else:
            #raise Exception('Reference checksums do not match disk')
        return

    def parse(self, limit=None):
        """

        :param limit:
        :return:
        """
        if limit is not None:
            logger.info("Only parsing first %s rows fo each file", str(limit))

        logger.info("Parsing files...")

        if self.testOnly:
            self.testMode = True
        self.label_hash = {'gene': {}, 'disease': {}}
        self._process_pathways(limit)
        self._process_diseases(limit)
        self._process_genes(limit)
        self._process_kegg_disease2gene(limit)
        #TODO: Finish omim2gene
        #self._process_omim2gene(limit)
        self._process_omim2disease(limit)
        self._process_genes_kegg2ncbi(limit)
        self._process_ortholog_classes(limit)

        for f in ['hsa_orthologs', 'mmu_orthologs', 'rno_orthologs','dme_orthologs','dre_orthologs','cel_orthologs']:
            file = '/'.join((self.rawdir, self.files[f]['file']))
            self._process_orthologs(file, limit)



        logger.info("Finished parsing")

        self.load_bindings()

        logger.info("Found %d nodes", len(self.graph))
        return

    def _process_pathways(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing pathways")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['pathway']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (pathway_id, pathway_name) = row

                if self.testMode and pathway_id not in self.test_ids['pathway']:
                    continue

                pathway_id = 'KEGG:'+pathway_id.strip()
                # Add the pathway as a class.
                gu.addClassToGraph(g, pathway_id, pathway_name)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with pathways")
        return

    def _process_diseases(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing diseases")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (disease_id, disease_name) = row

                if self.testMode and disease_id not in self.test_ids['disease']:
                    continue

                disease_id = 'KEGG:'+disease_id.strip()
                # Add the disease as a class.
                gu.addClassToGraph(g, disease_id, disease_name)
                if disease_id not in self.label_hash['disease']:
                    self.label_hash['disease'][disease_id] = disease_name

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with diseases")
        return

    def _process_genes(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing genes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['hsa_genes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, gene_name) = row

                if self.testMode and gene_id not in self.test_ids['genes']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                # Add the gene as a class.
                gu.addClassToGraph(g, gene_id, gene_name)
                if gene_id not in self.label_hash['gene']:
                    self.label_hash['gene'][gene_id] = gene_name

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with genes")
        return

    def _process_ortholog_classes(self, limit=None):
        """

        :param limit:
        :return:
        """

        logger.info("Processing ortholog classes")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['ortholog_classes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (orthology_class_id, orthology_class_name) = row

                if self.testMode and orthology_class_id not in self.test_ids['ortholog_classes']:
                    continue

                #FIXME: What's the proper route for this?
                # The orthology class is essentially a KEGG gene ID that is species agnostic.
                # Add the ID and label as a class. Would it be considered a gene as well?
                orthology_class_id = 'KEGG:'+orthology_class_id.strip()
                orthology_symbols = re.sub(';.*','',orthology_class_name)
                orthology_description = re.sub('.*;','',orthology_class_name)
                #FIXME: Problem if there is more than one symbol?
                #FIXME: Is there a designated type for these orthology classes?
                gu.addClassToGraph(g, orthology_class_id, orthology_symbols, None, orthology_description)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with ortholog classes")
        return


    def _process_orthologs(self, raw, limit=None):
        """
        This method maps orthologs for a species to the KEGG orthology classes.

        Triples created:
        <gene_id> is a class
        <orthology_class_id> is a class

        <assoc_id> has subject <gene_id>
        <assoc_id> has object <orthology_class_id>
        :param limit:
        :return:
        """

        logger.info("Processing orthologs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        gu = GraphUtils(curie_map.get())
        #raw = ('/').join((self.rawdir, self.files['ortholog_classes']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, orthology_class_id) = row

                #if self.testMode and orthology_id not in self.test_ids['ortholog_classes']:
                    #continue

                orthology_class_id = 'KEGG:'+orthology_class_id.strip()
                gene_id = 'KEGG:'+gene_id.strip()

                # note that the panther_id references a group of orthologs, and is not 1:1 with the rest
                assoc_id = self.make_id(''.join((gene_id, orthology_class_id)))
                #ortho = OrthologyAssoc()
                rel = OrthologyAssoc.ortho_rel['orthologous']
                # add the association and relevant nodes to graph
                assoc = OrthologyAssoc(assoc_id, gene_id, orthology_class_id, None, None)
                assoc.setRelationship(rel)
                assoc.loadAllProperties(g)  # FIXME inefficient

                # add gene and orthology class to graph; assume labels will be taken care of elsewhere
                gu.addClassToGraph(g, gene_id, None)
                gu.addClassToGraph(g, orthology_class_id, None)
                assoc.addAssociationToGraph(g)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with orthologs")
        return

    def _process_kegg_disease2gene(self, limit=None):
        """
        This method creates an association between diseases and their associated genes.

        Triples created:
        <alternate_locus> is an Individual
        <alternate_locus> has type <variant_locus>
        <alternate_locus> is an allele of  <gene_id>

        <assoc_id> has subject <disease_id>
        <assoc_id> has object <gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing disease to gene")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['disease_gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, disease_id) = row

                if self.testMode and gene_id not in self.test_ids['']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                disease_id = 'KEGG:'+disease_id.strip()

                # Make an association ID.
                assoc_id = self.make_id((disease_id+gene_id))

                # we actually want the association between the gene and the disease to be via an alternate locus
                # not the "wildtype" gene itself.
                # so we make an anonymous alternate locus, and put that in the association.
                #FIXME: Should we use the self.makeID here,
                # or is it better to see the gene & disease IDs for this alt_locus?
                alt_locus = '_'+gene_id+'-'+disease_id+'VL'
                alt_label = self.label_hash['gene'].get(gene_id)
                disease_label = self.label_hash['disease'].get(disease_id)
                if alt_label is not None and alt_label != '':
                    alt_label = 'some variant of '+alt_label+' that causes '+disease_label
                else:
                    alt_label = None
                gu.addIndividualToGraph(g, alt_locus, alt_label, geno.genoparts['variant_locus'])
                geno.addAlleleOfGene(alt_locus, gene_id)
                # Add the disease to gene relationship.
                assoc = G2PAssoc(assoc_id, alt_locus, disease_id, None, None)
                assoc.loadAllProperties(g)
                assoc.addAssociationToGraph(g)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with disease to gene")
        return

    def _process_omim2gene(self, limit=None):
        """
        This method creates an association between omim diseases and their associated genes.

        Triples created:
        <alternate_locus> is an Individual
        <alternate_locus> has type <variant_locus>
        <alternate_locus> is an allele of  <gene_id>

        <assoc_id> has subject <disease_id>
        <assoc_id> has object <gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing OMIM to gene")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['omim2gene']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (gene_id, disease_id) = row

                if self.testMode and gene_id not in self.test_ids['']:
                    continue

                gene_id = 'KEGG:'+gene_id.strip()
                disease_id = 'KEGG:'+disease_id.strip()

                # Make an association ID.
                assoc_id = self.make_id((disease_id+gene_id))

                # we actually want the association between the gene and the disease to be via an alternate locus
                # not the "wildtype" gene itself.
                # so we make an anonymous alternate locus, and put that in the association.
                #FIXME: Should we use the self.makeID here,
                # or is it better to see the gene & disease IDs for this alt_locus?
                alt_locus = '_'+gene_id+'-'+disease_id+'VL'
                alt_label = self.label_hash['gene'].get(gene_id)
                disease_label = self.label_hash['disease'].get(disease_id)
                if alt_label is not None and alt_label != '':
                    alt_label = 'some variant of '+alt_label+' that causes '+disease_label
                else:
                    alt_label = None
                gu.addIndividualToGraph(g, alt_locus, alt_label, geno.genoparts['variant_locus'])
                geno.addAlleleOfGene(alt_locus, gene_id)
                # Add the disease to gene relationship.
                assoc = G2PAssoc(assoc_id, alt_locus, disease_id, None, None)
                assoc.loadAllProperties(g)
                assoc.addAssociationToGraph(g)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with disease to gene")
        return

    def _process_omim2disease(self, limit=None):
        """
        This method will map the KEGG disease IDs to the corresponding OMIM disease IDs.
        Currently this only maps KEGG diseases and OMIM diseases that have a 1:1 mapping.
        Triples created:
        <kegg_disease_id> is a class
        <omim_disease_id> is a class
        <kegg_disease_id> hasEquivalentClass <omim_disease_id>
        :param limit:
        :return:
        """

        logger.info("Processing KEGG disease to OMIM disease mappings.")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        omim_disease_hash = {}
        kegg_disease_hash = {}
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['omim2disease']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (omim_disease_id, kegg_disease_id, link_type) = row

                if self.testMode and omim_disease_id not in self.test_ids['']:
                    continue

                kegg_disease_id = 'KEGG:'+kegg_disease_id.strip()
                omim_disease_id = re.sub('omim','OMIM',omim_disease_id)

                if omim_disease_id not in omim_disease_hash:
                    omim_disease_hash[omim_disease_id] = [kegg_disease_id]
                else:
                    omim_disease_hash[omim_disease_id].append(kegg_disease_id)
                    #print('Found additional kegg disease id for '+omim_disease_id)

                if kegg_disease_id not in kegg_disease_hash:
                    kegg_disease_hash[kegg_disease_id] = [omim_disease_id]
                else:
                    kegg_disease_hash[kegg_disease_id].append(omim_disease_id)
                    #print('Found additional omim disease id for '+kegg_disease_id)

                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        #Now process the disease hashes and only process 1:1 omim disease:KEGG disease entries.
        for omim_disease_id in omim_disease_hash:
            if len(omim_disease_hash[omim_disease_id]) == 1:
                kegg_disease_id = ('').join(omim_disease_hash.get(omim_disease_id))
                if len(kegg_disease_hash[kegg_disease_id]) == 1:
                    #print(kegg_disease_id+'_'+omim_disease_id)
                    gu.addClassToGraph(g, kegg_disease_id, None)
                    gu.addClassToGraph(g, omim_disease_id, None)
                    gu.addEquivalentClass(g, kegg_disease_id, omim_disease_id)

        logger.info("Done with KEGG disease to OMIM disease mappings.")
        return

    def _process_genes_kegg2ncbi(self, limit=None):
        """
        This method will map the KEGG human gene IDs to the corresponding NCBI Gene IDs.
        Triples created:
        <kegg_gene_id> is a class
        <ncbi_gene_id> is a class
        <kegg_gene_id> equivalentClass <ncbi_gene_id>
        :param limit:
        :return:
        """

        logger.info("Processing KEGG gene IDs to NCBI gene IDs")
        if self.testMode:
            g = self.testgraph
        else:
            g = self.graph
        line_counter = 0
        geno = Genotype(g)
        gu = GraphUtils(curie_map.get())
        raw = ('/').join((self.rawdir, self.files['ncbi']['file']))
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (kegg_gene_id, ncbi_gene_id) = row

                if self.testMode and kegg_gene_id not in self.test_ids['']:
                    continue

                # Adjust the NCBI gene ID prefix.
                ncbi_gene_id = re.sub('ncbi-geneid','NCBIGene',ncbi_gene_id)
                kegg_gene_id = 'KEGG:'+kegg_gene_id

                # Adding the KEGG gene ID to the graph here is redundant, unless there happens to be
                # additional gene IDs in this table not present in the genes table.
                gu.addClassToGraph(g, kegg_gene_id, None)
                gu.addClassToGraph(g, ncbi_gene_id, None)
                gu.addEquivalentClass(g, kegg_gene_id, ncbi_gene_id)


                if (not self.testMode) and (limit is not None and line_counter > limit):
                    break

        logger.info("Done with KEGG gene IDs to NCBI gene IDs")
        return