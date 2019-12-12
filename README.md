[![Build Status](https://travis-ci.org/smsaladi/interpro-simple-query.svg?branch=master)](https://travis-ci.org/smsaladi/interpro-simple-query)

Simple Interpro Query
=====================

This is a simple python implementation of interproscan's precalculated
query for proteins with precalculated domains. Should work for proteins 
that are a part of the current release of Uniprot. This might be easier to
use than interproscan itself if you know your sequences are all in Uniprot.

It is an alternative to downloading the complete domain matches
[`match_complete.xml.gz`](https://www.ebi.ac.uk/interpro/download/)
and parsing this file for proteins of interest.

Try:

```shell
pip install -r requirements.txt
python query_interpro_precalc.py --max_per_query 25 --sleep 10 proteins.fasta > proteins.ipr.ndjson &
```

Beware: The script can break at any point if Interpro changes the structure
of their API.

See code for implementation and other details. There are links to relevant bits
of interproscan to better understand Interpro's API pattern (see docstring).
