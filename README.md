# OncoHIT
Framework for the estimation of Intra-Tumor Heterogeneity

To estimate Intra-Tumor Heterogeneity (ITH) we use two methods:

- EXPANDS (Andor et al., 2014)

- PhyloWGS (Deshwar et al., 2015)

The inputs to EXPANDS are somatic single nucleotide variants (SNVs) and copy number estimates. 

To obtain somatic SNVs we use Mutect (Cibulskis et al., 2013), following an approach simiar to the Broad best practices.

For details on how we run this analysis, you can see a separate [vcall repository](https://github.com/comicsfct/vcall).
