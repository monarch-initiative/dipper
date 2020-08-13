"""
A fork of rdflib rdf2dot utility, see
https://rdflib.readthedocs.io/en/stable/_modules/rdflib/tools/rdf2dot.html#rdf2dot

We apply the formatliteral function to labels, but otherwise the code is the same

This is necessary for variants with HGVS primary labels that contain characters
that need to be url encoded (<, >)

Also replaces cgi.escape with html.escape

TO DO make a PR
"""

import rdflib
import rdflib.extras.cmdlineutils

import html
import collections

from rdflib import XSD

LABEL_PROPERTIES = [rdflib.RDFS.label,
                    rdflib.URIRef("http://purl.org/dc/terms/title"),
                    rdflib.URIRef("http://xmlns.com/foaf/0.1/name"),
                    rdflib.URIRef("http://www.w3.org/2006/vcard/ns#fn"),
                    rdflib.URIRef("http://www.w3.org/2006/vcard/ns#org")
                    ]

XSDTERMS = [
    XSD[x] for x in (
        "anyURI", "base64Binary", "boolean", "byte", "date",
        "dateTime", "decimal", "double", "duration", "float", "gDay", "gMonth",
        "gMonthDay", "gYear", "gYearMonth", "hexBinary", "ID", "IDREF",
        "IDREFS", "int", "integer", "language", "long", "Name", "NCName",
        "negativeInteger", "NMTOKEN", "NMTOKENS", "nonNegativeInteger",
        "nonPositiveInteger", "normalizedString", "positiveInteger", "QName",
        "short", "string", "time", "token", "unsignedByte", "unsignedInt",
        "unsignedLong", "unsignedShort")]

EDGECOLOR = "blue"
NODECOLOR = "black"
ISACOLOR = "black"


def rdf2dot(g, stream, graph_opts=None):
    """
    Convert the RDF graph to DOT
    writes the dot output to the stream
    """

    fields = collections.defaultdict(set)
    nodes = {}

    if graph_opts is None: graph_opts = {}

    def node(x):

        if x not in nodes:
            nodes[x] = "node%d" % len(nodes)
        return nodes[x]

    def label(x, g):

        for labelProp in LABEL_PROPERTIES:
            l = g.value(x, labelProp)
            if l:
                return formatliteral(l, g)

        try:
            return g.namespace_manager.compute_qname(x)[2]
        except:
            return x

    def formatliteral(l, g):
        v = html.escape(l)
        if l.datatype:
            return '&quot;%s&quot;^^%s' % (v, qname(l.datatype, g))
        elif l.language:
            return '&quot;%s&quot;@%s' % (v, l.language)
        return '&quot;%s&quot;' % v

    def qname(x, g):
        try:
            q = g.compute_qname(x)
            return q[0] + ":" + q[2]
        except:
            return x

    def color(p):
        return "BLACK"

    stream.write("digraph { \n node [ fontname=\"DejaVu Sans\" ] ; \n")

    for option, value in graph_opts.items():
        stream.write(" {}={} \n".format(option, value))

    for s, p, o in g:
        sn = node(s)
        if p == rdflib.RDFS.label:
            continue
        if isinstance(o, (rdflib.URIRef, rdflib.BNode)):
            on = node(o)
            opstr = "\t%s -> %s [ color=%s, label=< <font point-size='10' " + \
                    "color='#336633'>%s</font> > ] ;\n"
            stream.write(opstr % (sn, on, color(p), qname(p, g)))
        else:
            fields[sn].add((qname(p, g), formatliteral(o, g)))

    for u, n in list(nodes.items()):
        stream.write("# %s %s\n" % (u, n))
        f = ["<tr><td align='left'>%s</td><td align='left'>%s</td></tr>" %
             x for x in sorted(fields[n])]
        opstr = "%s [ shape=none, color=%s label=< <table color='#666666'" + \
                " cellborder='0' cellspacing='0' border='1'><tr>" + \
                "<td colspan='2' bgcolor='grey'><B>%s</B></td></tr><tr>" + \
                "<td href='%s' bgcolor='#eeeeee' colspan='2'>" + \
                "<font point-size='10' color='#6666ff'>%s</font></td>" + \
                "</tr>%s</table> > ] \n"
        stream.write(opstr % (n, NODECOLOR, label(u, g), u, u, "".join(f)))

    stream.write("}\n")
