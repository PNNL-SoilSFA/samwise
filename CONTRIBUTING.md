# Contributing to SAMWISE

Thank you for your interest in improving SAMWISE. Contributions to workflows,
helper scripts, documentation, examples, and reproducibility practices are
welcome.

Please read this guide before opening an issue or pull request. Contributions
must also follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Before you start

- Search existing issues and pull requests before opening a new one.
- For a substantial change, open an enhancement issue first so the proposed
  direction can be discussed.
- Do not include human or environmental metadata, private sequencing data,
  credentials, API keys, tokens, or generated analysis results in a commit.
- Do not report a suspected vulnerability in a public issue. Follow
  [SECURITY.md](SECURITY.md) instead.

## Fork-based workflow

External contributors should work from a personal fork. Direct pushes to the
upstream repository are reserved for maintainers.

```bash
git clone https://github.com/<your-account>/samwise.git
cd samwise
git remote add upstream https://github.com/PNNL-SoilSFA/samwise.git
git switch -c <type>/<short-description>
```

Keep the branch focused on one logical change. Update it from `upstream/main`
when needed, then push the branch to your fork and open a pull request against
`PNNL-SoilSFA/samwise:main`.

Maintainers may use a different local workflow, but should preserve the same
review and security expectations.

## Commits

Use [Conventional Commits](https://www.conventionalcommits.org/) for commit
messages. Common examples include:

```text
docs: clarify Module 2 inputs
feat: add grouped assembly option
fix: handle missing manifest rows
security: avoid exposing agent credentials
```

Keep commits small and explain the reason for non-obvious implementation
choices. Pull requests from external contributors are normally squash-merged
into `main`.

## Pull requests

A good pull request should:

- Explain what changed and why.
- Link an issue when one exists. Issue links are preferred but are not required
  for small changes.
- Identify affected workflows, scripts, documentation, or outputs.
- Include relevant commands, datasets, versions, and environment details used
  for validation.
- Update documentation and examples when behavior or interfaces change.
- Avoid unrelated formatting or generated-output changes.

Every pull request requires at least one maintainer review. It must also pass
the repository's required GitHub Advanced Security checks before approval and
merge. Do not disable, bypass, or conceal a security check failure.

The project is still defining its complete automated testing and validation
policy. Until that policy is finalized, describe the validation you performed
and clearly identify anything you could not run.

## Workflow and data practices

- Run workflows against representative, non-sensitive test data.
- Keep large databases, workflow results, temporary directories, and local
  environments outside the repository unless a specific small fixture is
  intentionally part of an example.
- Record relevant tool and database versions when reporting scientific results
  or workflow behavior.
- Treat downloaded third-party databases and bundled reference files according
  to their own licenses and terms.
- Use local environment files for credentials. Never commit `.env` files or
  replace example placeholders with real secrets.

## Documentation contributions

Documentation is built with Zensical. From the repository root, install
Zensical in a local virtual environment and run:

```bash
zensical serve
```

Open `http://localhost:8000` to preview changes. The published site is built by
GitHub Actions from the `main` branch.

## Questions

For questions that are not security-sensitive, open an issue or discussion in
the repository. For suspected vulnerabilities, use the private process in
[SECURITY.md](SECURITY.md).
