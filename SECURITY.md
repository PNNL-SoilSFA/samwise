# Security Policy

## Reporting a vulnerability

Please do not report suspected vulnerabilities in a public GitHub issue,
discussion, or pull request.

Send a private report to `csoc@pnnl.gov` with:

- A clear description of the issue and its potential impact
- The affected file, workflow, version, or commit
- Reproduction steps or a minimal proof of concept, if safe to provide
- Any relevant logs, screenshots, or suggested mitigation
- Your preferred method of contact and whether you would like attribution

Please avoid including private sequencing data, credentials, API keys, or other
sensitive material in email. If sensitive evidence is required, ask for a
secure transfer method.

## Response process

The security team will acknowledge reports, investigate and coordinate a fix,
and determine disclosure timing. Response and remediation timelines may depend
on severity, reproducibility, affected dependencies, and organizational
requirements.

## Supported versions

The latest released version and the `main` development branch are the primary
targets for security fixes. Older releases may require an upgrade before a fix
is available.

## Security expectations for contributors

- Never commit credentials, tokens, private data, or generated secrets.
- Do not add code that weakens validation or bypasses security controls without
  documenting and reviewing the reason.
- Use pinned or versioned dependencies where practical.
- Report accidentally exposed secrets privately and rotate them immediately.
