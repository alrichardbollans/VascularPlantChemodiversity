## A Dataset Of Phytochemistry and Species Chemodiversity

This dataset uses `phytochempy` [1] to download all compound occurrences for vascular plants in WikiData [2] and KNApSAcK [3]. 
Using the compound occurrences, the chemodiversity of species is calculated using a range of measures from [1]. 
The finalised lists of compound occurrences and diversity are given in `outputs`.

To cite this dataset, please use the associated Zenodo archive:

` `



## Technical Details

All names are resolved to the World Checklist of Vascular Plants v14 using `wcvpy` [4].

For the final datasets, records with no associated chemical ID (InChIKey, SMILES or CAS ID) and records with unresolved plant names are excluded.
 Duplicates of chemical occurrences are removed, identified based on standardised SMILES keys and accepted species.

The measurements of diversity are sensitive to sampling effort, please read [1] for more in depth detail on this. 
The rarification process from [1] is extremely computationally expensive for this size dataset and is not included here.
 As seen in [1], it is likely that rarified FAD, MFAD and APWD resemble unrarified APWD and that rarified H, Hbc and G resemble unrarified versions.


### References & Acknowledgements

[1] Richard-Bollans, Adam, Eliot Jan-Smith, Daniele Silvestro, and Melanie-Jayne R. Howes. ‘Beyond Phylogeny: Phytochemical Diversity as a Unique Metric for Biodiversity in the Gentianales’. New Phytologist 248, no. 6 (2025). https://doi.org/10.1111/nph.70653.

[2] Denny Vrandečić and Markus Krötzsch, ‘Wikidata: A Free Collaborative Knowledgebase’, Communications of the ACM 57, no. 10 (2014): 78–85.

[3] Farit Mochamad Afendi, Taketo Okada,, Mami Yamazaki, Aki-Hirai-Morita, Yukiko Nakamura,
Kensuke Nakamura, Shun Ikeda, Hiroki Takahashi, Md. Altaf-Ul-Amin, Latifah, Darusman, Kazuki
Saito, Shigehiko Kanaya, “KNApSAcK Family Databases: Integrated Metabolite-Plant Species
Databases for Multifaceted Plant Research,” Plant Cell Physiol., 53, e1(1-12), (2012). doi:
10.1093/pcp/pcr165.

[4] Richard-Bollans, Adam. ‘wcvpy’. Zenodo, 2026. https://doi.org/10.5281/zenodo.14774384.

The developers acknowledge Research Computing at the James Hutton Institute for providing computational resources and technical support for the 'UK’s
Crop Diversity Bioinformatics HPC' (BBSRC grants BB/S019669/1 and BB/X019683/1), use of which has contributed to the development of the model used in
this analysis.
