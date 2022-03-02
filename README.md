# OncoHIT
Framework for the estimation of Intra-Tumor Heterogeneity

To estimate Intra-Tumor Heterogeneity (ITH) we use two methods:

- EXPANDS (Andor et al., 2014)

- PhyloWGS (Deshwar et al., 2015)

The inputs to EXPANDS are somatic single nucleotide variants (SNVs) and copy number estimates. 

To obtain somatic SNVs we use Mutect (Cibulskis et al., 2013), following an approach similar to the Broad best practices.

For details on how we obtain somatic SNVs you can see the separate [vcall repository](https://github.com/comicsfct/vcall).

To obtain copy number estimates we use cnvkit (Talevich et al., 2016). Namely, we use the cnvkit batch command using alignments from tumor samples (obtained from the variant calling procedure) against a pool of alignments of the normal samples.

