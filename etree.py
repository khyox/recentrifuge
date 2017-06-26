#from lxml import etree
from xml.dom import minidom
import xml.etree.ElementTree as etree

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = etree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")

krona = etree.Element('krona', attrib={'collapse':'false', 'key':'true'})
ET = etree.ElementTree(krona)

attributes = etree.SubElement(krona, 'attributes', {'magnitude':'grams'})
attribute1 = etree.SubElement(attributes, 'attribute',
                              {'display':'Grams'})
attribute1.text = 'grams'
attribute2 = etree.SubElement(attributes, 'attribute',
                              {'display':'% Daily Value'})
attribute2.text = 'dva'
dss = etree.SubElement(krona, 'datasets')
ds1 = etree.SubElement(dss, 'dataset')
ds1.text = 'Brand X'
ds2 = etree.SubElement(dss, 'dataset')
ds2.text = 'Brand Y'
ds3 = etree.SubElement(dss, 'dataset')
ds3.text = 'Brand Z'

root = etree.SubElement(krona, 'node', name='Granola serving')
gr = etree.SubElement(root, 'grams')
val = etree.SubElement(gr, 'val')
val.text = '55'
val = etree.SubElement(gr, 'val')
val.text = '55'
val = etree.SubElement(gr, 'val')
val.text = '85'

prot = etree.SubElement(root, 'node', name='Protein', href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=7711')
gr = etree.SubElement(prot, 'grams')
val = etree.SubElement(gr, 'val')
val.text = '5'
val = etree.SubElement(gr, 'val')
val.text = '6'
val = etree.SubElement(gr, 'val')
val.text = '8'

test = etree.SubElement(root, 'node', name='Test')
gr = etree.SubElement(test, 'grams')
val = etree.SubElement(gr, 'val')
val.text = '25'
val = etree.SubElement(gr, 'val')
val.text = '0'
val = etree.SubElement(gr, 'val')
val.text = '16'
gr = etree.SubElement(test, 'dva')
val = etree.SubElement(gr, 'val')
val.text = '125'
val = etree.SubElement(gr, 'val')
val.text = '0'
val = etree.SubElement(gr, 'val')
val.text = '416'

subtest = etree.SubElement(test, 'node', name='subtest')
gr = etree.SubElement(subtest, 'grams')
val = etree.SubElement(gr, 'val')
val.text = '5'
val = etree.SubElement(gr, 'val')
val.text = '0'
val = etree.SubElement(gr, 'val')
val.text = '6'


#xml.etree.ElementTree.SubElement(parent, tag, attrib={}, **extra)

print(prettify(krona))

#ET.write('xml_test.xml')
with open('xml_test.xml', 'w') as xml_file:
    xml_file.write(prettify(krona))
