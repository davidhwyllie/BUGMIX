# TBMIX
Study mixtures of TB sequences in VCF files

This includes four main classes:
* vcfStore- parses VCF files, computing minor allele frequencies and storing them in on-disc highly indexed sqlite databases, one per vcf file.
* mixtureReader - extracts information from one vcfStore'd item.  Various methods are available extracting useful information 
* multiMixtureReader - extracts information from multiple vcfStore'd items.  

* summaryStore - a helper class used by the mixtureReader and multiMixtureReader classes to persist summary information generated to disc.  
             - this ensures that computations, which may be expensive, only occur once.
             - the data is accessible from the RDBMS level but the persistence is transparent if using the mixture & multiMixtureReader classes.
             
Uses cases addressed concern the identification of mixed samples of TB.
