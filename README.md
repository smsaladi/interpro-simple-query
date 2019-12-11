Simple Interpro Query
=====================

This is a simple python implementation to query interpro for
proteins with precalculated domains. Should work for proteins 
that are a part of the current release of Uniprot.

This is an alternative to downloading the complete domain matches
`match_complete.xml` and parsing this file for proteins of interest.

Can break at any point if Interpro changes the structure of their API.
Based on the interproscan precalculated query code (see comments) 

