import sys
import xml.parsers.expat

do_print = False

def start_element(name, attrs):
    global do_print
    if name == 'Well' and attrs['row'] == 'A' and attrs['col'] == '2':
        do_print = True
    if do_print:
        print name, attrs

def end_element(name):
    global do_print
    if name == 'Well':
        do_print = False

for f in sys.argv[1:]:
    try:
        print f
        parser = xml.parsers.expat.ParserCreate()
        parser.returns_unicode = False
        parser.StartElementHandler = start_element
        parser.EndElementHandler = end_element
        parser.ParseFile(open(f))
    except Exception, e:
        pass
