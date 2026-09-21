# Optional AI agent

SAMWISE includes a beta OpenAI-based AI agent for interrogating genomes and
their metabolisms from files in an output directory. It is separate from the
Nextflow workflows, is not activated by default, and is optional.

For the agent's original setup notes, see the
[`agent/README.md`](https://github.com/PNNL-SoilSFA/samwise/blob/main/agent/README.md).

From the `agent/` directory, create the environment:

```bash
bash setupAgentOpenai.sh
conda activate langchain-chat-openai
```

Copy the example environment file and set the API key and user settings
locally:

```bash
cp .env.example .env
```

Run the agent from `agent/` with:

```bash
python chatOpenai.py
```

The agent is an analysis assistant, not a Nextflow orchestrator: it does not
start or manage SAMWISE workflow runs. It can inspect files under the roots in
`ALLOWED_ROOTS` and can write an explicitly requested, approved module only to
the locations configured in `MODULE_WRITE_ROOTS`; those capabilities are
disabled or restricted by local configuration. It is usable only after the
environment, provider settings, and Python launcher have been explicitly
configured and activated. SAMWISE's core workflows do not require the agent.

`.env` is local configuration and must not be committed. Never commit API
keys, access tokens, credentials, private datasets, or generated conversation
memory. Follow the [contribution guide](https://github.com/PNNL-SoilSFA/samwise/blob/main/CONTRIBUTING.md)
and [security policy](https://github.com/PNNL-SoilSFA/samwise/blob/main/SECURITY.md)
for repository security practices.
