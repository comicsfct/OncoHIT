# OncoHIT
Framework for the estimation of Intra-Tumor Heterogeneity

To estimate Intra-Tumor Heterogeneity (ITH) we use two methods:

- EXPANDS (Andor et al., 2014)

- PhyloWGS (Deshwar et al., 2015)

To infer clones and their frequencies with EXPANDS we run expands_clones.R

The inputs to EXPANDS are somatic single nucleotide variants (SNVs) and copy number estimates. 

To obtain somatic SNVs we use Mutect (Cibulskis et al., 2013), following an approach similar to the Broad best practices.

For details on how we obtain somatic SNVs you can see the separate [vcall repository](https://github.com/comicsfct/vcall).

To obtain copy number estimates we use cnvkit (Talevich et al., 2016). Namely, we use the cnvkit batch command using alignments from tumor samples (obtained from the variant calling procedure) against a pool of alignments of the normal samples.

For PhyloWGS we also integrate purity and tumor cell fraction information from Absolute (Carter et al., 2015). 


As an extra, we also include scripts to calculate other measures of heterogeneity. This includes very simple measures of tumor mutation burden using SNVs or CNVs. We also provide a script to estimate the Tumor Heterogeneity index TH (Oh et al., 2018) and the MATH score (Mroz er al., 2015). 

Finally, we also include a script to process somatic variants and copy number variants from one (or several samples) to be used with pyclone-vi (Gillis and Roth, 2020).
