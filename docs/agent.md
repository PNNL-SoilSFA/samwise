# Optional AI agent

SAMWISE includes a beta AI agent for interrogating files written to an output
directory. The agent is separate from the Nextflow workflows and is not
activated by default.

To configure it, follow the setup instructions in `agent/readme.txt` and use a
local `.env` file containing your own provider configuration. Never commit API
keys, access tokens, credentials, private datasets, or generated conversation
memory to the repository.

The agent is optional. SAMWISE's core workflows can be used without installing
or configuring it.
