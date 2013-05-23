DNA Electrostatics in R
============================================================

Fast algorithm of calculation of DNA electrostatics was developed in 1999. Since then we have looked through regulatory elements of many genomes and show that electrostatical properties of regulatory DNA sites could explain a lot of features of DNA-protein recognition process.

All analysis in that papers was done with Matlab, which is not free. Active development of R and Bioconductor project especially make us to migrate our code to R/Bioconductor, so it will be much more accessible to bioinformatics community.

INSTALLATION
-----------------------------------------------------------
The **reldna** library is pure  R code so it should work on any platform where R is installed.
Installation process is really straightforward:

1. Download [latest version](https://github.com/promodel/reldna/raw/master/distr/reldna_0.0-9.20130208.tar.gz) of the library
2. Use standard R installation procedure
   * Chechk that all dependencies are installed in your system
   * Install library 
```
R CMD INSTALL reldna_THE_VERSION_YOU_VE_DOWNLOADED.tar.gz
```

REFERENCES
-----------------------------------------------------------
There is long list of publications about electrostatic properties of regulatory DNA and bacterial promoters in particular:

1. KAMZOLOVA, SG; BESKARAVAINY, PM; SOROKIN, AA. Electrostatic Properties of Promoter and Nonpromoter Sites in T7 Bacteriophage Genome. JOURNAL OF BIOMOLECULAR STRUCTURE & DYNAMICS, 2009,  26:6 p. 881-882.
2. Kamzolova, S. G., Sorokin, A. A., Osipov, A. A., Beskaravainy, P. M. Electrostatic Map of Bacteriophage T7 Genome. Comparative Analysis of Electrostatic Properties of sigma(70)-Specific T7 DNA Promoters Recognized by RNA-Polymerase of Escherichia coli. Biofizika, 2009, 54:6 p. 975-983.
3. Sorokin, A.A., et al., Analysis of the distribution of the nucleotide sequence and electrostatic potential of the Escherichia coli genome. Biofizika, 2007. 52(2): p. 223-227.
4. Sorokin, A.A., et al., Electrostatic properties of promoters recognized by E. coli RNA polymerase Esigma-70. Journal of bioinformatics and computational biology, 2006. 4(2): p. 455-67.
5. Sorokin, A.A., et al. Electrostatic Map of E. coli Genome DNA. Specific Features of Electrostatic Potential of Promoter and Nonpromoter Regions. in The 14th Albany Conversation. 2005. University of Albany, Albany, NY.
6. Sorokin, A.A., et al. Oligonucleotide Analysis of E. coli Promoters Recognized by σ70-RNA Polymerase. in The 14th Albany Conversation. 2005. University of Albany, Albany, NY.
7. Kamzolova, S.G., et al., [Some principles in the organization of σ70-specific promoters on the E. coli genome on the basis of electrostatic patterns of promoter DNA]. Biofizika, 2005. 50(3): p. 444-9.
8. Kamzolova, S.G., et al., Electrostatic potentials of E.coli genome DNA. J Biomol Struct Dyn, 2005. 23(3): p. 341-5.
9. Kamzolova, S.G., et al., Comparative analysis of electrostatic patterns for promoter and non-promoter DNA in E. coli, in Bioinformatics of Genome Regulation And Structure II, N. Kolchanov, et al., Editors. 2005, Springer. p. 67-74.
10. Sorokin, A.A., et al. Electrostatic Properties of Promoter DNA. Approach to Classification Analysis. in 13th Albany Conversation. 2003. University of Albany, Albany, NY.
11. Kamzolova, S.G., et al. Specific Features of Electrostatic Potential of E.coli Ribosomal Promoters. in 13th Albany Conversation. 2003. University of Albany, Albany, NY.
12. Sorokin, A.A., et al. The quest for new forms of promoter determinants. Relationship of promoter nucleotide sequences to their electrostatic potential distribution. in 12th Albany Conversation. 2001. University of Albany, Albany, NY.
13. Kamzolova, S.G., et al., RNA polymerase--promoter recognition. Specific features of electrostatic potential of "early" T4 phage DNA promoters. J Biomol Struct Dyn, 2000. 18(3): p. 325-34.
14. Polozov, R.V., et al., Electrostatic potentials of DNA. Comparative analysis of promoter and nonpromoter nucleotide sequences. J Biomol Struct Dyn, 1999. 16(6): p. 1135-43.
