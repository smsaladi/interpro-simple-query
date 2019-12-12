"""

This is a simple python implementation to query interpro for proteins with
precalculated domains. Should work for proteins that are a part of the current
release of Uniprot. If a query sequence has no precalculated entry, then it
will not be present in the output. If there are duplicate sequences sent in
the same query, only one of the two will be returned.

No sequence cleanup is performed, e.g. an `*` at the end of a protein,
say for a stop codon, will result in nothing being found, since
the Uniprot entry wouldn't include this.


Output is given on standard out as newline-delimited JSON. The sequence itself
is added to each entry, for convenience.

API reference gleaned from from interproscan codebase:

* set Sequence
https://github.com/ebi-pf-team/interproscan/blob/99af1b48/core/model/src/main/java/uk/ac/ebi/interpro/scan/model/Protein.java#L333-L354

* calculate MD5
https://github.com/ebi-pf-team/interproscan/blob/99af1b48/core/model/src/main/java/uk/ac/ebi/interpro/scan/model/Protein.java#L534-L546

* query precalculated
https://github.com/ebi-pf-team/interproscan/blob/99af1b48/core/precalcmatches/precalc-match-client/src/main/java/uk/ac/ebi/interpro/scan/precalc/client/MatchHttpClient.java#L154-L165


Author: Shyam Saladi (saladi@caltech.edu)
Date: December 2019
License: MIT

"""

import argparse
import hashlib
from itertools import chain, islice
from collections import OrderedDict
from urllib3.util.retry import Retry

import requests
from requests.adapters import HTTPAdapter
import Bio.SeqIO
import xmltodict
import ndjson

INTERPRO_PRECALC_URL = 'https://www.ebi.ac.uk/interpro/match-lookup/matches/'


def chunk(iterable, n):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))

def md5seq(seq):
    """Calculates the md5 of a sequence. No error checking
    """
    seq = seq.upper()
    return hashlib.md5(seq.encode('utf-8')).hexdigest().upper()


# Set up requests with retries and backoff
# https://stackoverflow.com/a/35504626/2320823
sess = requests.Session()
retries = Retry(total=5,
                backoff_factor=0.1,
                status_forcelist=[ 500, 502, 503, 504 ])
sess.mount('http://', HTTPAdapter(max_retries=retries))

def query_interpro(seqs):
    """Makes a single request for one or more sequences
    """
    if not isinstance(seqs, list):
        seqs = [seqs]
    
    seqs_md5 = OrderedDict([(md5seq(s), s) for s in seqs])

    payload = [('md5', m) for m, s in seqs_md5.items()]
    response = sess.post(INTERPRO_PRECALC_URL, data=payload)

    data = xmltodict.parse(response.content.decode('utf-8'))
    data = data['kvSequenceEntryXML']['matches']['match']

    # Only one hit was returned
    if isinstance(data, dict):
        data = [data]

    # Add protein sequence to response
    for i, elem in enumerate(zip(data, seqs_md5)):
        ipr, seq_md5 = elem
        data[i]['seq'] = seqs_md5[ipr['proteinMD5']]

    return data

def cast_to_seq(seqgen):
    """Convert generator of Bio.SeqRecords to sequence strings
    """
    for rec in seqgen:
        yield str(rec.seq)
    return

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_fn')
    parser.add_argument('--format', default='fasta')
    parser.add_argument('--max_per_query', default=10, type=int)

    args = parser.parse_args()

    gen = cast_to_seq(Bio.SeqIO.parse(args.seq_fn, args.format))

    for batch in chunk(gen, args.max_per_query):
        data = query_interpro(batch)
        # print output as NDjson
        print(ndjson.dumps(data))

    return


# These shouldn't fall out of Uniprot (their domain/HMM annotations may evolve)
ecoli_secy = \
    "MAKQPGLDFQSAKGGLGELKRRLLFVIGALIVFRIGSFIPIPGIDAAVLAKLLEQQRGTIIEMFNMFSGGALSRASIFALGIMPYISASIIIQLLTVVHPTLAEIKKEGESGRRKISQYTRYGTLVLAIFQSIGIATGLPNMPGMQGLVINPGFAFYFTAVVSLVTGTMFLMWLGEQITERGIGNGISIIIFAGIVAGLPPAIAHTIEQARQGDLHFLVLLLVAVLVFAVTFFVVFVERGQRRIVVNYAKRQQGRRVYAAQSTHLPLKVNMAGVIPAIFASSIILFPATIASWFGGGTGWNWLTTISLYLQPGQPLYVLLYASAIIFFCFFYTALVFNPRETADNLKKSGAFVPGIRPGEQTAKYIDKVMTRLTLVGALYITFICLIPEFMRDAMKVPFYFGGTSLLIVVVVIMDFMAQVQTLMMSSQYESALKKANLKGYGR"

bsubt_secy = \
    "MFKTISNFMRVSDIRNKIIFTLLMLIVFRIGAFIPVPYVNAEALQAQSQMGVFDLLNTFGGGALYQFSIFAMGITPYITASIIIQLLQMDVVPKFTEWSKQGEVGRRKLAQFTRYFTIVLGFIQALGMSYGFNNLANGMLIEKSGVSTYLIIALVLTGGTAFLMWLGEQITSHGVGNGISIIIFAGIVSSIPKTIGQIYETQFVGSNDQLFIHIVKVALLVIAILAVIVGVIFIQQAVRKIAIQYAKGTGRSPAGGGQSTHLPLKVNPAGVIPVIFAVAFLITPRTIASFFGTNDVTKWIQNNFDNTHPVGMAIYVALIIAFTYFYAFVQVNPEQMADNLKKQGGYIPGVRPGKMTQDRITSILYRLTFVGSIFLAVISILPIFFIQFAGLPQSAQIGGTSLLIVVGVALETMKQLESQLVKRNYRGFMKN"

def test_md5seq():
    assert md5seq(ecoli_secy) == '4CD214BD2268CC5AB7F4398B9CBC6CE0'
    assert md5seq(ecoli_secy) == md5seq(ecoli_secy.lower())
    assert md5seq(ecoli_secy) != md5seq(ecoli_secy + 'X')
    return

def test_query_interpro():
    # Makes sure queries resolve (can't do much other than that, since HMMs change...)
    # and checks that the output format is consistent/as expected

    # Check single sequence query
    data = query_interpro(ecoli_secy)
    assert type(data) == list
    assert len(data[0]['hit']) > 0
    assert 'seq' in data[0].keys()

    # Check multiple, distinct sequence query
    seqs = [ecoli_secy, bsubt_secy]
    data = query_interpro(seqs)
    assert len(data) == 2
    assert type(data) == list
    for s, d in zip(seqs, data):
        assert len(d['hit']) > 0
        assert d['seq'] == s
        assert d['proteinMD5'] == md5seq(s)

    # If two identical sequences are sent, you only get one response
    data = query_interpro([ecoli_secy, ecoli_secy])
    assert len(data) == 1
    assert type(data) == list

    # With an 'X', the first sequence won't resolve (not precalculated)
    # only the proper one is returned with data
    data = query_interpro([ecoli_secy + 'X', ecoli_secy])
    assert len(data) == 1
    assert len(data[0]['hit']) > 0
    assert data[0]['seq'] == ecoli_secy
    assert data[0]['proteinMD5'] == md5seq(ecoli_secy)

    return


if __name__ == '__main__':
    main()
