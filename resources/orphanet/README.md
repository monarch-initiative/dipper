
## Orphanet XML

To survey:

    xmlstarlet el -a  ../../raw/orphanet/en_product6.xml | sort -u > orphanet6.xpath

    xpath2dot.awk -v ORIENT="UD"  en_product6.xpath > en_product6.gv

    dot -Tpng en_product6.gv > en_product6.png

    ![schema](./en_product6.png)


note: 2018 Fall, there is no `GeneList` element in the XML
although it is refered to in the code and exists in the dipper/test/resources/orphanet
examples.

    grep GeneList en_product6.xpath

goose-egg

