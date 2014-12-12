from rdflib import BNode, ConjunctiveGraph, Literal, RDF, Graph
from rdflib.namespace import FOAF, DC, RDFS, OWL

import urllib, csv, os, time
from urllib import request
from datetime import datetime
from stat import *
import hashlib


core_bindings = {'dc': DC, 'foaf': FOAF, 'rdfs': RDFS}



class Source:
    def __init__(self,args=[]):
        self.name = args
        self.path = ""
        self.graph = ConjunctiveGraph()
        self.triple_count = 0
        self.outdir = 'out'
        self.rawdir = 'raw'
        return

    # here's where we will have a configured
    #dictionary of prefix to uri bindings for the graph
    def load_core_bindings(self):
        self.graph.bind("dc", DC)
        self.graph.bind("foaf", FOAF)
        self.graph.bind("rdfs", RDFS)
        self.graph.bind('owl', OWL)

        return

    def load_bindings(self):
        #should each source provide it's own binding dictionary?
        #or should it be a universal one, like the current context file?
        return

    def fetch(self):
        return

    def parse(self):

        #self.addnode("abc:123","testnode",'type:12345','test description here')
        return

    def write(self):
        for s, p, o in self.graph:
            print(s, p, o)
        self.graph.serialize(format="turtle")
        return

    def whoami(self):
        print("I am ",self.name)
        return

    def addnode(self, id, label, type, description=None):
        mynode = BNode()
        #sself.graph.add((mynode,DC['id'],id))
        self.graph.add((mynode, RDFS['label'], Literal(label)))
        if description is not None:
            self.graph.add((mynode, DC['description'], Literal(description)))
        #self.graph.add((mynode, RDF.type, type))

        return

    def make_id(self,long_string):
        #for now, this will do md5.  probably not the best long-term solution
        #note others available: md5(), sha1(), sha224(), sha256(), sha384(), and sha512()
        byte_string = long_string.encode("utf-8")
        return hashlib.md5(byte_string).hexdigest()


    def checkIfRemoteIsNewer(self,remote,local):
        print("Checking if remote file is newer...")
        #check if local file exists
        #if no local file, then remote is newer
        if (not os.path.exists(local)):
            print("File does not exist locally")
            return True
        #get remote file details
        d = urllib.request.urlopen(remote)
        last_modified = d.info()['Last-Modified']
        print("Remote file date:",last_modified)
        #"2008-11-10 17:53:59"
        #dt_obj = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
        #Thu, 07 Aug 2008 16:20:19 GMT
        dt_obj = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")

        size=d.info()['Content-Length']

        #get local file details
        st = os.stat(local)
        print("Local file date:",datetime.utcfromtimestamp(st[ST_CTIME]))

        #check date on local vs remote file
        if (dt_obj > datetime.utcfromtimestamp(st[ST_CTIME])):
            #check if file size is different
            if (st[ST_SIZE] != size):
                print("Newer file exists on remote server")
                return True
            else:
                print("Remote file has same filesize--will not download")
        return False

    def fetch_from_url(self,remotefile,localfile):
        if (self.checkIfRemoteIsNewer(remotefile, localfile)):
            print("Fetching from ", remotefile)
            # TODO url verification, etc
            annotation_file = urllib.request
            annotation_file.urlretrieve(remotefile, localfile)
            print("Finished.  Wrote file to", localfile)
        else:
            print("Using existing file.")

        st = os.stat(localfile)
        print("file size:", st[ST_SIZE])
        print("file created:", time.asctime(time.localtime(st[ST_CTIME])))
        return

#store.add((donna, RDF.type, FOAF.Person))
#store.add((donna, FOAF.nick, Literal("donna", lang="foo")))
#store.add((donna, FOAF.name, Literal("Donna Fales")))

