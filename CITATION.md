# Citation

Please cite the SAMWISE release used in your work. The preferred citation will
be finalized after the project team confirms the author list, affiliations, and
ORCID identifiers.

## Software citation

> **SAMWISE: [confirmed project title].** [Confirmed authors]. Version
> `[VERSION]`. Zenodo. DOI: `[ZENODO_DOI]`.

Replace the placeholders after the corresponding Zenodo release is available.
The DOI should point to the version-specific Zenodo record rather than only to
the repository landing page.

## Metadata to confirm

- Project title: `[CONFIRM TITLE]`
- Authors in citation order: `[CONFIRM AUTHORS]`
- Institutional affiliations: `[CONFIRM AFFILIATIONS]`
- ORCID identifiers: `[CONFIRM ORCIDs]`
- Cited version: `[VERSION]`
- Release date: `[RELEASE DATE]`
- Zenodo DOI: `[ZENODO_DOI]`

## BibTeX template

```bibtex
@software{samwise_[YEAR]_[REPOSITORY],
  author  = {[CONFIRM AUTHORS]},
  title   = {{SAMWISE: [CONFIRM PROJECT TITLE]}},
  version = {[VERSION]},
  year    = {[YEAR]},
  publisher = {Zenodo},
  doi     = {[ZENODO_DOI]},
  url     = {https://github.com/PNNL-SoilSFA/samwise}
}
```

Do not replace the DOI placeholder until Zenodo has assigned the official DOI.
The repository may later add a `CITATION.cff` file after the team confirms the
same metadata.

## Citing results and dependencies

Software citation does not replace citation of the scientific methods,
reference databases, and third-party tools used in an analysis. Cite those
dependencies according to their own documentation and licenses, and record
their versions alongside the SAMWISE release.
