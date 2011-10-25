import sys
import xml.parsers.expat

print_all = 0
depth = 0

do_print = False

def start_element(name, attrs):
    global do_print
    if name == 'Filters':
        do_print = True
    if do_print:
        print name, attrs

def end_element(name):
    if name == 'Filters':
        raise ValueError('foo')

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
