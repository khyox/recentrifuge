"""
Classes and functions directly related with Krona.

"""
# pylint: disable=not-an-iterable
import csv
import html
import os
import subprocess
import xml.etree.ElementTree as ETree
from xml.dom import minidom
from typing import List, Dict, NewType, Any, Optional

from recentrifuge.config import Filename, Sample, Scoring, Chart
from recentrifuge.config import JSLIB, HTML_SUFFIX
from recentrifuge.config import yellow, red
from recentrifuge.stats import SampleStats

# Type annotations
# pylint: disable=invalid-name
Attrib = NewType('Attrib', str)  # Refers to Krona attributes not XML ones
Elm = ETree.Element
# pylint: enable=invalid-name

# Predefined constants
COUNT = Attrib('count')
UNASSIGNED = Attrib('unassigned')
TID = Attrib('tid')
RANK = Attrib('rank')
SCORE = Attrib('score')
HREFBASE_TAX = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id='
HREFBASE_GEN = 'http://amigo.geneontology.org/amigo/term/'
LOGO_RCF = '/img/logo-rcf-mini.uri'
LOGO_RGF = '/img/logo-rgf-mini.uri'
HIDDEN = '/img/hidden.uri'
LOADING = '/img/loading.uri'
FAVICON = '/img/favicon.uri'

# Define encoding dialect for TSV files expected by Krona
csv.register_dialect('krona', csv.get_dialect('unix'), delimiter='\t',
                     quoting=csv.QUOTE_NONE)


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
        """Wrapper for creating a meaningful Krona node.

        For details, please consult:
        https://github.com/marbl/Krona/wiki/Krona-2.0-XML-Specification
        """
        subnode = self.sub(parent, 'node',
                           {'name': name,
                            'href': f'https://www.google.es/search?q={name}'})
        count_node = self.sub(subnode, COUNT)
        counts: Dict[Sample, str] = {sample: values[COUNT][sample]
                                     for sample in self.samples}
        for sample in self.samples:
            counts_value: Optional[str] = counts[sample]
            # Save space (warning! Empty tag instead of 0 inside <val></val>)
            if counts_value is not None and int(counts_value) == 0:
                counts_value = None
            self.sub(count_node, 'val', None, counts_value)
        if values.get(UNASSIGNED) and any(values[UNASSIGNED].values()):
            # Avoid including and save space if all the unassigned values are 0
            unassigned_node = self.sub(subnode, UNASSIGNED)
            unassigned: Dict[Sample, str] = {sample: values[UNASSIGNED][sample]
                                             for sample in self.samples}
            for sample in self.samples:
                unassigned_value: Optional[str] = unassigned[sample]
                # Save space (warning! Empty and not 0 after <val>)
                if unassigned_value is not None and int(unassigned_value) == 0:
                    unassigned_value = None
                self.sub(unassigned_node, 'val', None, unassigned_value)
        if values.get(TID):
            tid_node = self.sub(subnode, TID)
            self.sub(tid_node, 'val',
                     {'href': values[TID]},
                     values[TID])
        if values.get(RANK):
            rank_node = self.sub(subnode, RANK)
            self.sub(rank_node, 'val', None, values[RANK])
        if values.get(SCORE):
            score_node = self.sub(subnode, SCORE)
            scores: Dict[Sample, str] = {sample: values[SCORE][sample]
                                         for sample in self.samples}
            for sample in self.samples:
                self.sub(score_node, 'val', None, scores[sample])
        return subnode

    @staticmethod
    def to_pretty_string(element: Elm):
        """Return a pretty-printed XML string for the Element."""
        raw_string = ETree.tostring(element,
                                    encoding='unicode',
                                    method='xml',
                                    short_empty_elements=False,
                                    )
        re_parsed = minidom.parseString(raw_string)
        pretty = re_parsed.toprettyxml(indent='  ')
        pretty = html.unescape(pretty)
        return pretty.split('\n', 1)[-1]  # Remove the XML 1.0 tag

    def __init__(self,
                 samples: List[Sample],
                 num_raw_samples: int = None,
                 stats: Dict[Sample, SampleStats] = None,
                 min_score: float = 0.0,
                 max_score: float = 1.0,
                 scoring: Scoring = Scoring.SHEL,
                 chart: Chart = Chart.TAXOMIC,
                 ) -> None:
        """
        Args:
            samples: List of samples in the set
            num_raw_samples: Number of raw samples (not from cross-analysis)
            min_score: minimum expected score
            max_score: maximum expected score
            chart: Type of chart (taxonomic, genomic...)
        """
        # Type declaration
        self.krona: Elm
        self.krona_tree: ETree.ElementTree
        self.attributes: Elm
        self.samples: List[Sample]
        self.datasets: Elm

        self.chart: Chart = chart

        # Dummy dict if stats not provided
        if stats is None:
            stats = {}

        # Select for taxonomic or genomic chart
        iden: str
        hrefbase: str
        display: str
        if self.chart == Chart.TAXOMIC:
            iden = 'TaxID'
            hrefbase = HREFBASE_TAX
            if scoring is Scoring.SHEL:
                display = 'Score (avg)'
            elif scoring is Scoring.LENGTH:
                display = 'Read length (avg)'
            elif scoring is Scoring.LOGLENGTH:
                display = 'Length (log10)'
            elif scoring is Scoring.NORMA:
                display = 'Score/Length (%)'
            elif scoring is Scoring.LMAT:
                display = 'LMAT score (avg)'
            elif scoring is Scoring.CLARK_C:
                display = 'CLARK conf (%)'
            elif scoring is Scoring.CLARK_G:
                display = 'CLARK gamma (avg)'
            elif scoring is Scoring.KRAKEN:
                display = 'Kmer coverage (%)'
            elif scoring is Scoring.GENERIC:
                display = 'Generic score (avg)'
            else:
                raise Exception(red('\nERROR!'),
                                f'Unknown Scoring "{scoring}"')
        elif self.chart == Chart.GENOMIC:
            iden = 'GenID'
            hrefbase = HREFBASE_GEN
            display = 'Score from Eval'
        else:
            raise Exception(f'ERROR! Unknown Chart "{self.chart}"')

        # Set root of KronaTree
        self.krona = ETree.Element('krona',  # type: ignore
                                   attrib={'collapse': 'true',
                                           'key': 'true',
                                           'chart': str(chart)})

        # Set attributes
        self.attributes = ETree.SubElement(self.krona, 'attributes',
                                           {'magnitude': 'count'})
        # # Set Count attribute
        self.sub(self.attributes, 'attribute',
                 {'display': 'Count', 'dataAll': 'members',
                  'tip': 'Number of reads assigned to this and child taxa'},
                 'count')
        # # Set Unassigned attribute
        self.sub(self.attributes, 'attribute',
                 {'display': 'Unassigned', 'dataNode': 'members',
                  'tip': 'Number of reads assigned specifically to this taxon'},
                 'unassigned')
        # # Set Id attribute
        self.sub(self.attributes, 'attribute',
                 {'display': iden, 'mono': 'true', 'hrefBase': hrefbase,
                  'tip': 'Taxonomic identifier'},
                 'tid')
        # # Set Rank attribute
        self.sub(self.attributes, 'attribute',
                 {'display': 'Rank', 'mono': 'true',
                  'tip': 'Taxonomic rank/level'},
                 'rank')
        # # Set confidence/score attribute
        self.sub(self.attributes, 'attribute',
                 {'display': display,
                  'tip': 'Averaged score of reads assigned'
                         ' to this and child taxa'},
                 'score')

        # Set datasets
        self.samples = samples
        self.datasets = ETree.SubElement(self.krona, 'datasets',
                                         {'rawSamples': f'{num_raw_samples}'})
        for sample in self.samples:
            if sample in stats:
                self.sub(self.datasets, 'dataset',
                         stats[sample].to_krona(), sample)
            else:
                self.sub(self.datasets, 'dataset', {}, sample)

        # Set color
        self.color = self.sub(self.krona, 'color',
                              {'attribute': 'score',
                               'hueStart': '0',
                               'hueEnd': '300',
                               'valueStart': f'{min_score:.1f}',
                               'valueEnd': f'{max_score:.1f}',
                               'default': 'true'},
                              ' ')  # Krona: Avoid empty-element tag

        super(KronaTree, self).__init__(self.krona)

    def __repr__(self):
        return self.to_pretty_string(self.krona)

    def tofile(self,
               filename: Filename,
               pretty: bool = False,
               ) -> None:
        """
        Write KronaTree in 'plain' or 'pretty' XML.

        Args:
            filename: the name of the XML output file.
            pretty: this parameter controls the layout of the XML code
                so that it is human readable for True (use for debug
                only because it uses a lot more of space and also has
                empty tags which are currently not supported by Krona)
                and machine readable for False (default, saves space).

        Returns: None

        """
        with open(filename, 'w') as xml_file:
            if pretty:
                xml_file.write(self.to_pretty_string(self.krona))
            else:
                self.write(xml_file,
                           encoding='unicode',
                           xml_declaration=False,
                           method='xml',
                           short_empty_elements=False,
                           )

    def tohtml(self,
               filename: Filename,
               pretty: bool = False,
               ) -> None:
        """
        Write Krona HTML.

        Args:
            filename: the name of the HTML output file.
            pretty: this parameter controls the layout of the XML code
                so that it is human readable for True (use for debug
                only because it uses a lot more of space and also has
                empty tags which are currently NOT SUPPORTED BY KRONA)
                and machine readable for False (default, saves space).

        Returns: None

        """

        # Warn about use of pretty option
        if pretty:
            print(yellow(f'\nWARNING! Pretty XML uses empty tags which are'
                         f' UNSUPPORTED by Krona-JS!'))
            print(yellow(f'WARNING! Prepare for unexpected HTML results!'))

        # Read aux files
        path = os.path.dirname(os.path.realpath(__file__))
        with open(path + HIDDEN, 'r') as file:
            hidden_image = file.read()
        with open(path + LOADING, 'r') as file:
            loading_image = file.read()
        with open(path + FAVICON, 'r') as file:
            favicon = file.read()
        path_logo: str
        if self.chart == Chart.TAXOMIC:
            path_logo = path + LOGO_RCF
        elif self.chart == Chart.GENOMIC:
            path_logo = path + LOGO_RGF
        else:
            raise Exception(f'ERROR! Unknown Chart "{self.chart}"')

        with open(path_logo, 'r') as file:
            logo = file.read()
        with open(f'{path}/{JSLIB}', 'r') as file:
            script = file.read()

        # Set root of HTML doc
        html_root = ETree.Element(  # type: ignore
            'html', attrib={'xmlns': 'http://www.w3.org/1999/xhtml',
                            'xml:lang': 'en',
                            'lang': 'en'})
        # Prepare HTML file
        head = self.sub(html_root, 'head')
        self.sub(head, 'meta', {'charset': 'utf-8'})
        self.sub(head, 'link', {'rel': 'shortcut icon',
                                'href': favicon})
        self.sub(head, 'link', {'rel': 'stylesheet',
                                'href': 'https://fonts.googleapis.com/css?family=Ubuntu'})
        self.sub(head, 'script', {'id': 'notfound'},
                 'window.onload=function(){document.body.innerHTML=""}')
        self.sub(head, 'script',
                 {'language': 'javascript', 'type': 'text/javascript'},
                 script)  # Include javascript
        body = self.sub(html_root, 'body')
        self.sub(body, 'img', {'id': 'hiddenImage',
                               'src': hidden_image,
                               'style': 'display:none'})
        self.sub(body, 'img', {'id': 'loadingImage',
                               'src': loading_image,
                               'style': 'display:none'})
        self.sub(body, 'img', {'id': 'logo',
                               'src': logo,
                               'style': 'display:none'})
        self.sub(body, 'noscript', None,
                 'Javascript must be enabled to view this page.')

        div = self.sub(body, 'div', {'style': 'display:none'})
        div.append(self.krona)  # Include specific XML from samples
        # Write the HTML file
        with open(filename, 'w') as html_file:
            html_file.write(
                '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n')  # pylint: disable=line-too-long
            if pretty:
                html_file.write(self.to_pretty_string(html_root))
            else:
                html_file.write(ETree.tostring(html_root,
                                               encoding='unicode',
                                               method='html',
                                               short_empty_elements=False,
                                               )
                                )


def krona_from_xml(xmlfile: Filename,
                   htmlfile: Filename = Filename('Output' + HTML_SUFFIX),
                   ):
    """Generate the Krona html file calling ktImportXML."""
    subprc = ["ktImportXML"]
    subprc.append(xmlfile)
    subprc.extend(["-o", htmlfile])
    try:
        subprocess.run(subprc, check=True)
    except subprocess.CalledProcessError:
        print('\n\033[91mERROR!\033[0m ktImportXML: ' +
              'returned a non-zero exit status (Krona plot built failed)')
