import pandas as pd
from dipper.graph.RDFGraph import RDFGraph
from dipper.models.Model import Model

columns = ['variant', 'variant_label', 'variant_type',
           'phenotype','relation', 'source', 'evidence', 'dbxref']

data =  [
     ['ClinVarVariant:254143', 'C326F', 'SO:0000694',
      'HP:0000504','RO:0002200', 'PMID:12503095', 'ECO:0000220',
      'dbSNP:886037891']
]

# Initialize graph and model
graph = RDFGraph()
model = Model(graph)

# Read file
dataframe = pd.DataFrame(data=data, columns=columns)

for index, row in dataframe.iterrows():
   # Add the triple ClinVarVariant:254143 RO:0002200 HP:0000504
   # RO:0002200 is the has_phenotype relation
   # HP:0000748 is the phenotype 'Inappropriate laughter'
   model.addTriple(row['variant'], row['relation'], row['phenotype'])

   # The addLabel method adds a label using the rdfs:label relation
   model.addLabel(row['variant'], row['variant_label'])

   # addType makes the variant an individual of a class,
   # in this case SO:0000694 'SNP'
   model.addType(row['variant'], row['variant_type'])

   # addXref uses the relation OIO:hasDbXref
   model.addXref(row['variant'], row['dbxref'])

   # Serialize the graph as turtle
   print(graph.serialize(format='turtle').decode("utf-8"))
