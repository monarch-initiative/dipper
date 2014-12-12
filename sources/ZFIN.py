import csv
from utils import pysed
import os, datetime
from datetime import datetime
from stat import *


from sources.Source import Source
from models.Assoc import Assoc
from models.Genotype import Genotype
from models.Dataset import Dataset
from rdflib import Namespace


class ZFIN(Source):
    GENOTYPES_URL = "http://zfin.org/downloads/genotype_features.txt"
    GENOTYPES_FILE = "genotype_features.txt"

    namespaces = {
        'ZP': 'http://purl.obolibrary.org/obo/ZP_',
        'ZFIN': 'http://zfin.org/',
        'GENO': 'http://purl.obolibrary.org/obo/GENO_',
    }

    def __init__(self):
        Source.__init__(self, 'zfin')
        self.outfile = self.outdir + '/' + self.name + ".ttl"
        self.datasetfile = self.outdir + '/' + self.name + '_dataset.ttl'
        self.genofile = ('/').join((self.rawdir, self.GENOTYPES_FILE))
        print("Setting outfile to", self.outfile)
        self.namespaces.update(Assoc.namespaces)

        self.dataset = Dataset('zfin', 'ZFIN', 'http://www.zfin.org')
        self.dataset.setFileAccessUrl(self.GENOTYPES_URL)
        #zfin versions are set by the date of download
        st = os.stat(self.genofile)

        #self.dataset.setVersion(datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d"),datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d_%H:%M:%S"))
        self.dataset.setVersion(datetime.utcfromtimestamp(st[ST_CTIME]).strftime("%Y-%m-%d"))

        return

    def load_bindings(self):
        self.load_core_bindings()
        for k in self.namespaces.keys():
            v=self.namespaces[k]
            self.graph.bind(k, Namespace(v))

        return

    def fetch(self):
        self.fetch_from_url(self.GENOTYPES_URL, self.genofile)

        self.scrub()

        return

    def scrub(self):
        # scrub file of the oddities where there are "\" instead of empty strings
        pysed.replace("\\\\", '', self.genofile)

        return

    # here we're reading and building a full named graph of this resource, then dumping it all at the end
    # we can investigate doing this line-by-line later
    # supply a limit if you want to test out parsing the head X lines of the file
    def parse(self, limit=None):
        Source.parse(self)
        if (limit is not None):
            print("Only parsing first", limit, "rows")
        line_counter = 0

        self._process_genotype_features(self.genofile, self.outfile, self.graph, limit)
        self.load_bindings()

        #self._process_phenotype_tab(self.rawfile,self.outfile,self.g,limit)


        filewriter = open(self.outfile, 'w')
        print(self.graph.serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        filewriter = open(self.datasetfile,'w')
        print(self.dataset.getGraph().serialize(format="turtle").decode(), file=filewriter)
        filewriter.close()

        print("Wrote", self.triple_count, "things")
        return

    def _process_genotype_features(self, raw, out, g, limit=None):
        line_counter = 0
        with open(raw, 'r', encoding="utf8") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1

                (genotype_id, genotype_name, genotype_unique_name, allele_id, allele_name, allele_ab,
                 allele_type, allele_disp_type, gene_symbol, gene_id, zygosity,
                 construct_name, construct_id, other) = row

                genotype_id = 'http://www.zfin.org/' + genotype_id.strip()
                geno = Genotype(genotype_id, genotype_name)

                #reassign the allele_type to a proper GENO or SO class
                allele_type = self._map_allele_type_to_geno(allele_type)

                allele_id = 'http://www.zfin.org/' + allele_id.strip()
                geno.addAllele(allele_id, allele_name, allele_type)

                if (gene_id is not None and gene_id.strip() != ''):
                    gene_id = 'http://www.zfin.org/' + gene_id.strip()
                    geno.addGene(gene_id, gene_symbol)

                    #if it's a transgenic construct, then we'll have to get the other bits
                    if (construct_id is not None and construct_id.strip() != ''):
                        construct_id = 'http://www.zfin.org/' + construct_id.strip()
                        geno.addAlleleDerivesFromConstruct(allele_id, construct_id)

                    #allele to gene
                    #FIXME i don't know why, but this shows up as ns1:GENO_0000440 instead of the URI
                    geno.addAlleleOfGene(allele_id, gene_id)


                #genotype has_part allele
                geno.addAlleleToGenotype(genotype_id, allele_id)
                #need to make some attributes of this relationship for allele zygosity within the genotype
                #or do we just make the allele_complement here, based on the zygosity?
                #allele_in_gene_id=self.make_id(genotype_id+allele_id+zygosity)
                #allele has_disposition zygosity?
                #                g.add(())

                if (limit is not None and line_counter > limit):
                    break

                self.graph = geno.getGraph().__iadd__(self.graph)
        return

    def _map_allele_type_to_geno(self, allele_type):
        type = None
        type_map = {
            'complex_substitution': 'SO:1000005',  #complex substitution
            'deficiency': 'SO:1000029',  #incomplete chromosome
            'deletion': 'SO:0000159',  #deletion
            'indel': 'SO:1000032',  #indel
            'insertion': 'SO:0000667',  #insertion
            'point_mutation': 'SO:1000008',  #point_mutation
            'sequence_variant': 'SO:0001060',  #sequence variant
            'transgenic_insertion': 'SO:0001218',  #transgenic insertion
            'transgenic_unspecified': 'SO:0000781',  #transgenic unspecified
            'transloc': 'SO:0000199',  #translocation
            #            'unspecified' : None
        }
        if (allele_type.strip() in type_map):
            type = 'http://purl.obolibrary.org/obo/' + type_map.get(allele_type)
#            print("Mapped: ", allele_type, "to", type)
        else:
            #TODO add logging
            print("ERROR: Allele Type (", allele_type, "not mapped")

        return type

