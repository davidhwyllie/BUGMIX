# BUGMIX
Allows storage and fast random access to read depths in VCF files derived from bacterial reference mapping.

BugMix includes four main classes:
* vcfStore- parses VCF files, computing minor allele frequencies and storing them in on-disc in indexed sqlite databases, one per vcf file.
* mixtureReader - extracts information from one vcfStore'd item.  Various methods are available extracting per-base, per-gene, or per arbitrary position set etc. 
* multiMixtureReader - extracts information from multiple vcfStore'd items.  

* summaryStore - a helper class used by the mixtureReader and multiMixtureReader classes to persist summary information generated in an RDBMS.
             - this ensures that computations, which may be expensive, only occur once.
             - the data is accessible from the RDBMS level but the persistence is transparent if using the mixture & multiMixtureReader classes.
             
Use cases addressed with this component concern 
- the identification of mixed samples of TB.
- identification of regions where mis-mapping may have occurred
