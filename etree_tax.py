"""KronaTree stuff"""

from typing import Dict, List
from typing import NewType, Any
from xml.dom import minidom
import xml.etree.ElementTree as ETree

# Type annotations
# pylint: disable=invalid-name
Filename = NewType('Filename', str)
Sample = NewType('Sample', str)
Attrib = NewType('Attrib', str)  # Refers to Krona attributes not XML ones
Elm = ETree.Element
# pylint: enable=invalid-name

# Predefined constants
COUNT = Attrib('count')
UNASSIGNED = Attrib('unassigned')
RANK = Attrib('rank')
SCORE = Attrib('score')


class KronaTree(ETree.ElementTree):
    """Kronified ElementTree."""

    @staticmethod
    def sub(parent: Elm,
            tag: str,
            attrib: Dict[str, str] = None,
            text: str = None,
            ) -> Elm:
        """Wrapper around ETree.SubElement."""
        if attrib:
            subelement = ETree.SubElement(parent, tag, attrib)
        else:
            subelement = ETree.SubElement(parent, tag)
        if text is not None:
            subelement.text = text
        return subelement

    def node(self,
             parent: Elm,
             name: str,
             values: Dict[Attrib, Any],
             ) -> Elm:
        """Wrapper for creating a meaningful Krona node."""
        subnode = self.sub(parent, 'node',
                           {'name': name,
                            'href': f'https://www.google.es/search?q={name}'})
        count_node = self.sub(subnode, COUNT)
        counts = {sample: values[COUNT][sample] for sample in self.samples}
        for sample in counts:
            self.sub(count_node, 'val', None, counts[sample])
        if values.get(RANK):
            rank_node = self.sub(subnode, RANK)
            self.sub(rank_node, 'val', None, values[RANK])
        if values.get(SCORE):
            score_node = self.sub(subnode, SCORE)
            scores: Dict[Sample, str] = {smpl: values[SCORE][smpl] for smpl in
                                         self.samples}
            for sample in scores:
                self.sub(score_node, 'val', None, scores[sample])
        return subnode

    @staticmethod
    def to_pretty_string(element: Elm):
        """Return a pretty-printed XML string for the Element."""
        raw_string = ETree.tostring(element, 'utf-8')
        re_parsed = minidom.parseString(raw_string)
        pretty = re_parsed.toprettyxml(indent='  ')
        return pretty.split('\n', 1)[-1]  # Remove the XML 1.0 tag

    def __init__(self,
                 samples: List[Sample],
                 ) -> None:
        # Type declaration
        self.krona: Elm
        self.krona_tree: ETree.ElementTree
        self.attributes: Elm
        self.samples: List[Sample]
        self.datasets: Elm

        # Set root of KronaTree
        self.krona = ETree.Element('krona',
                                   attrib={'collapse': 'false', 'key': 'true'})

        # Set attributes
        self.attributes = ETree.SubElement(self.krona, 'attributes',
                                           {'magnitude': 'count'})
        self.sub(self.attributes, 'attribute',
                 {'display': 'Count', 'dataAll': 'members'},
                 'count')
        self.sub(self.attributes, 'attribute',
                 {'display': 'Unassigned', 'dataNode': 'members'},
                 'unassigned')
        self.sub(self.attributes, 'attribute',
                 {'display': 'Rank', 'mono': 'true'},
                 'rank')
        self.sub(self.attributes, 'attribute',
                 {'display': 'Avg. confidence'},
                 'score')

        # Set datasets
        self.samples = samples
        self.datasets = ETree.SubElement(self.krona, 'datasets')
        for sample in self.samples:
            self.sub(self.datasets, 'dataset', {}, sample)

        # Set color
        self.color = self.sub(self.krona, 'color',
                              {'attribute': 'score',
                               'hueStart': '0',
                               'hueEnd': '120',
                               'valueStart': '0.0',
                               'valueEnd': '1.0',
                               'default': 'false'},
                              ' ')  # Krona: Avoid empty-element tag

        super(KronaTree, self).__init__(self.krona)

    def __repr__(self):
        return self.to_pretty_string(self.krona)

    def tofile(self,
              filename: Filename):
        """Write KronaTree in 'pretty' XML."""
        with open(filename, 'w') as xml_file:
            xml_file.write(self.to_pretty_string(self.krona))


def main():
    """Main entry point to code."""
    # Argument Parser Configuration

    s1 = Sample('S1')
    s2 = Sample('S2')

    samples: List[Sample] = [s1, s2, ]
    tree = KronaTree(samples)

    root = tree.node(tree.getroot(), 'root', {COUNT: {s1: '1000', s2: '2000'}})
    eukarya = tree.node(root, 'Eukaria',
                        {COUNT: {s1: '500', s2: '500'},
                         RANK: 'Superkingdom',
                         SCORE: {s1: '0.3', s2: '0.8'}})
    chordate = tree.node(eukarya, 'Chordate',
                         {COUNT: {s1: '250', s2: '125'},
                          RANK: 'Phylum',
                          SCORE: {s1: '0.2', s2: '0.9'}})

    print(tree)
    tree.tofile(Filename('xml_tax.xml'))


if __name__ == '__main__':
    main()
