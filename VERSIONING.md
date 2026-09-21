# Versioning and Releases

SAMWISE uses [Semantic Versioning](https://semver.org/) in the form
`MAJOR.MINOR.PATCH`. Release tags use the same version with a leading `v`, for
example `v1.2.0`.

## Version increments

### Major version

Increment the major version for breaking changes or a major rewrite. Breaking
changes include changes to workflow interfaces, input formats, output contracts,
or behavior that requires users to update existing commands or analyses.

Examples:

- Removing or renaming a workflow parameter
- Changing a manifest format incompatibly
- Replacing a workflow architecture in a way that changes its public behavior

### Minor version

Increment the minor version for backward-compatible feature revisions and
larger security patches. A minor release may add capabilities without requiring
existing users to change their commands.

Examples:

- Adding a new optional workflow parameter
- Adding a new auxiliary module
- Adding a larger security fix that does not break supported usage

### Patch version

Increment the patch version for small security updates and small changes that
do not change functionality. Backward-compatible changes may also use a patch
release when they are small, limited in scope, and amount to 500 lines of code
or fewer.

The line-count guideline is not a substitute for assessing behavior. A change
that affects compatibility or adds a meaningful feature should use the major or
minor category even if its implementation is short.

Examples:

- Correcting a documentation or packaging error
- Fixing a small bug without changing the supported interface
- Applying a small security hardening change

## Release process

Before creating a release, maintainers should:

1. Confirm that the version category matches the user-visible impact.
2. Review the changes since the previous release and update documentation.
3. Confirm that required security checks pass.
4. Create an annotated Git tag such as `v1.0.0` on the release commit.
5. Create a corresponding GitHub Release with notes for users.
6. Confirm that Zenodo has archived the GitHub Release and generated or updated
   the DOI record.
7. Update `CITATION.md` with the release version and DOI when the DOI is
   available.

Zenodo is connected to this repository and is expected to create the first DOI
when the first release is published. Until then, documentation should use the
placeholder shown in `CITATION.md` rather than inventing a DOI.

## Pre-releases

Pre-release versions may use SemVer identifiers such as `v1.0.0-rc.1` or
`v1.0.0-beta.1`. Pre-releases are not the default citation target unless a
release explicitly instructs users to cite one.

## Reproducibility

Release notes should identify the workflow and database versions that affect
results. Users should retain the SAMWISE version, commit or release tag, input
manifest, relevant command line, and external database versions with their
analysis records.
