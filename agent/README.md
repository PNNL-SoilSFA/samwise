# SAMWISE Query Agent - Setup and Usage Guide

This is a terminal chatbot that can answer questions and *do things* with files on disk: read documents, browse and search directories, analyze Python code, explain Nextflow pipelines, read Excel/R data files, design new Nextflow modules with your feedback, and remember your feedback across sessions. It's built with [LangChain](https://python.langchain.com/) and talks to a chat model (an LLM) through an OpenAI-compatible API endpoint — by default, PNNL's internal "AI Incubator" gateway.

This guide assumes no prior experience with Python packaging, conda, or LangChain. If a step doesn't apply to you (e.g. you're not on an HPC cluster), skip it.

---

## 1. What's actually in this repository

| File / folder | What it's for |
|---|---|
| `chatOpenai.py` | **The main program.** Run this to chat with the agent. Talks to any OpenAI-compatible API (PNNL AI Incubator, OpenAI, etc). |
| `chatOllama.py` + `ollama.sh` | An alternative that talks to a *locally hosted* model via [Ollama](https://ollama.com/), instead of a remote API. Most people should use `chatOpenai.py` instead — only use this if you specifically want a fully local/offline model. |
| `skills/` | The "tools" the agent can call — one Python file per capability (filesystem access, RAG document search, Excel, R `.rds` files, Nextflow pipeline analysis, Nextflow module design, feedback memory, etc). You don't need to run these directly. |
| `.env` | **Your personal configuration and credentials.** Not checked into git — you create this yourself (see below). |
| `.env.example` | A template for `.env` with every setting explained. Copy this to `.env` and fill in your values. |
| `requirementsOpenai.txt` | The exact list of Python packages `chatOpenai.py` and its skills need. |
| `setupAgentOpenai.sh` | One-time script that creates the conda environment and installs everything in `requirementsOpenai.txt`. |
| `run.sh` | Convenience script: activates the conda environment, double-checks packages are installed, and launches `chatOpenai.py`. |
| `docs/` | Default folder the document-search tool indexes. Empty until you put files in it. |
| `chroma_db/`, `module_designs/`, `feedback_memory.db`, `large_file_state.db` | Local storage the agent creates automatically the first time it runs. **You do not need to create these yourself** — see the FAQ below. |

---

## 2. One-time setup

### 2.1 Get the code and go to the project folder

```bash
cd /rcfs/projects/samwise/query_agent      # or wherever you cloned it
```

### 2.2 Create the conda environment and install dependencies

Run the setup script once:

```bash
bash setupAgentOpenai.sh
```

This creates a conda environment called `langchain-chat-openai` (Python 3.11) and installs every package listed in `requirementsOpenai.txt`: `langchain`, `langchain-openai`, `langchain-community`, `langchain-chroma`, `langchain-huggingface`, `langchain-text-splitters`, `chromadb`, `sentence-transformers`, `pypdf`, `docx2txt`, `python-dotenv`, `openpyxl`, `numpy`.

> **If you already created the environment before this fix and are
> missing packages** (`langchain-chroma`, `langchain-community`,
> `langchain-huggingface`, `openpyxl`, etc.), that's because an older
> version of this script only installed `langchain`, `langchain-openai`,
> and `python-dotenv`. Just re-run:
> ```bash
> conda activate langchain-chat-openai
> pip install -r requirementsOpenai.txt
> ```
> to pick up everything that was missing. `run.sh` now also does this
> automatically every time you launch the chat, so this shouldn't recur.

### 2.3 Create your `.env` file

The agent reads all of its configuration (API key, model name, which folders it's allowed to touch, etc.) from a file named `.env` in this same folder. It is **not** included in the repository (it contains secrets), so you must create it yourself:

```bash
cp .env.example .env
```

Then open `.env` in any text editor and fill in your values. At minimum you need:

```
BASE_URL=https://ai-incubator-api.pnnl.gov
MODEL_NAME=gpt-5.4-project
API_KEY=sk-your-real-key-here
```

Everything else in `.env.example` has a working default and is explained inline with comments — read through it once, it's short.

A few things that trip people up (see also the FAQ at the bottom):

* **Quotes around values are optional.** `ALLOWED_ROOTS=/rcfs/projects/samwise/` and `ALLOWED_ROOTS="/rcfs/projects/samwise/"` both work identically — don't worry about which one to use.
* **Don't type a double slash (`//`)** anywhere in a path, e.g. `FEEDBACK_DB_PATH=/some/path//feedback_memory.db` — that's just a typo, not special syntax. Use a single slash: `/some/path/feedback_memory.db`.
* **Folders and database files do not need to exist beforehand.** The agent creates any folder or `.db` file it needs automatically the first time it's used (see the FAQ below for details).

### 2.4 (Optional) Put documents in `docs/`

If you want to use the built-in document search tool, drop `.txt`, `.md`, `.pdf`, or `.docx` files into the `docs/` folder (or point `DOCS_DIR` in
`.env` at another folder). This is optional — everything else works without it.

---

## 3. Running the agent

Every time you want to chat with the agent:

```bash
bash run.sh
```

This activates the `langchain-chat-openai` conda environment, makes sure all required packages are installed, and starts `chatOpenai.py`.

Or, if you've already activated the environment yourself:

```bash
conda activate langchain-chat-openai
python chatOpenai.py
```

You'll see a startup banner, then a `You:` prompt. Type a question or request and press Enter. A few special commands:

| Command | Effect |
|---|---|
| `/reset` | Clear conversation history and start a fresh conversation. |
| `/exit` or `/quit` | Quit the chat. (`Ctrl+C` / `Ctrl+D` also work.) |

### Do you need an active SLURM session to run this?

**No**, not for `chatOpenai.py`. It only makes outbound HTTPS requests to the PNNL AI Incubator API (or whichever `BASE_URL` you configured) and, the first time it runs, downloads a small embedding model (`all-MiniLM-L6-v2`, a few hundred MB) from Hugging Face to run locally on CPU. Neither of those needs a GPU or a batch job — a login node or your own laptop is fine.

A SLURM/GPU session (and `ollama.sh`) is only relevant if you deliberately choose the **local Ollama** path (`chatOllama.py`) instead, e.g. because you want a fully offline model with no external API calls.


## 4. What the agent can do

The system prompt in `chatOpenai.py` tells the model which tool to reach for based on what you ask. In practice you can just ask naturally, e.g.:

* "What's in the `docs/` folder about X?" → searches indexed documents
* "List the files in /rcfs/projects/samwise/some_run" → browses that folder (must be under `ALLOWED_ROOTS`)
* "Explain what this Nextflow pipeline does: /path/to/pipeline" → summarizes `.nf` files and their processes
* "Read this .rds file for me" → loads an R data file
* "What's in this Excel workbook / find rows containing X" → reads/searches `.xlsx` files
* "This pipeline is missing a module for Y, can you design it?" → finds gaps in a Nextflow pipeline, drafts a module spec, and asks for your feedback before writing anything to disk
* "That last design isn't right, do X instead" → the agent records your feedback and revises the design; approved designs can later be written out as real `.nf` files with `write_module_to_file`

The agent also remembers feedback you give it (via `feedback_memory.db`) and will recall it in later conversations, even after restarting.

---

## 5. Frequently asked questions

### How do I find which `MODEL_NAME` / embedding model values are valid?

**`MODEL_NAME`** (the chat model) must be one of the model IDs the gateway in `BASE_URL` actually serves. For the PNNL AI Incubator gateway, you can list them yourself with your API key:

```bash
curl -s https://ai-incubator-api.pnnl.gov/v1/models \
  -H "Authorization: Bearer YOUR_API_KEY" | python -m json.tool
```

That returns a JSON list of model ids (e.g. `gpt-5.4-project`, `gpt-4o-project`, `claude-opus-4-5-20251101-v1-project`, `gemini-2.5-pro-project`, etc.) — use one of those for `MODEL_NAME`. You can also browse the same gateway's Swagger UI in a browser at `https://ai-incubator-api.pnnl.gov` and look under the "model management" section, or ask whoever manages the gateway for the current list.

**`EMBEDDING_MODEL` / `FEEDBACK_EMBED_MODEL`** are different — those are local [sentence-transformers](https://www.sbert.net/) models used for document/feedback search, downloaded from Hugging Face the first time they're used (no PNNL gateway involved). Any model name from the [sentence-transformers model list on Hugging Face](https://huggingface.co/models?library=sentence-transformers) works; the default, `all-MiniLM-L6-v2`, is small, fast, and good enough for most use cases. You generally don't need to change this.

### Do files and folders need to be created ahead of time, or does the agent make them on the fly?

**You do not need to pre-create anything.** Every folder and database file the agent uses is created automatically, on demand, the first time it's needed:

* `docs/`, `chroma_db/` (document search index)
* `module_designs/` (saved Nextflow module design drafts)
* `feedback_memory.db` and its parent folder (feedback memory database)
* `large_file_state.db` (bookkeeping for reading very large files)
* Any directory listed in `MODULE_WRITE_ROOTS` (where generated `.nf` files get written) — created the first time you write a module into it

The only thing you must create yourself is `.env` (from `.env.example`), and optionally add files to `docs/` if you want document search to have something to find.

### Do I need quotation marks around `ALLOWED_ROOTS` (or any other `.env` value)?

No — quotes are optional. `python-dotenv` (the library that reads `.env`) strips them either way, so `ALLOWED_ROOTS=/rcfs/projects/samwise/` and `ALLOWED_ROOTS="/rcfs/projects/samwise/"` are exactly equivalent. Use whichever is easier for you to read; it's harmless to leave the quotes in or take them out.

### Is `FEEDBACK_DB_PATH=/some/path//feedback_memory.db` (double slash) supposed to look like that?

No — that's a typo, not required syntax. A double slash still technically works on Linux (the OS collapses it), but it's confusing and easy to introduce by accident when copy-pasting; just use a single slash: `FEEDBACK_DB_PATH=/some/path/feedback_memory.db`.

### Do you need an active SLURM session? Should you?

No, and no — see "Do you need an active SLURM session to run this?" above. `chatOpenai.py` makes plain network calls; it doesn't need compute-node resources.

---

### I ran the script but it hangs or errors out talking to the model — what's wrong?

Most commonly this is one of:

1. **Missing packages.** If you see `ModuleNotFoundError: No module named 'langchain_chroma'` (or `langchain_community`, `langchain_huggingface`, `openpyxl`, etc.), your conda environment is missing packages that an older setup script didn't install. Fix it with:
   ```bash
   conda activate langchain-chat-openai
   pip install -r requirementsOpenai.txt
   ```
2. **Wrong or expired `API_KEY`.** Double check you copied the whole key (including the `sk-` prefix) with no extra spaces or line breaks, and that it's a key that's actually valid for the gateway in `BASE_URL`.
3. **Slow model + short timeout.** Some models (especially reasoning models such as `gpt-5.4-project`) can take well over a minute to respond, especially under load. `chatOpenai.py` now sets a 300-second default timeout (configurable via `REQUEST_TIMEOUT` in `.env`) to avoid this; if you still see `[Error contacting model: ... timed out...]`, try raising `REQUEST_TIMEOUT` further, or try a faster model like `gpt-4o-project` to confirm the rest of your setup is correct before switching back.
4. **Network/firewall.** The machine running the script needs outbound HTTPS access to `BASE_URL`. If you're behind a restrictive firewall or VPN, that request may never reach the gateway at all.

### Where do `ALLOWED_ROOTS` and `MODULE_WRITE_ROOTS` matter, and why are they separate?

`ALLOWED_ROOTS` controls which directories the agent is allowed to *read* from (list files, read file contents, index documents, read Excel/RDS files, run scripts). `MODULE_WRITE_ROOTS` separately controls which directories the agent is allowed to *write* generated Nextflow module files into. They're kept separate on purpose so that giving the agent read access to a large shared project directory doesn't automatically also give it permission to write files there — you opt into write access explicitly, and only for the specific output folders you choose.

---

## 6. Example walkthrough

A realistic session, assuming a SAMWISE-style pipeline output directory under `ALLOWED_ROOTS` (e.g. `/rcfs/projects/samwise/some_pipeline_run/`):

```
You: What Nextflow pipelines are under /rcfs/projects/samwise/some_pipeline_run?
  [calling tool: find_nextflow_pipelines({'path': '/rcfs/projects/samwise/some_pipeline_run'})]
Agent: Found main.nf, nextflow.config, and 6 modules under modules/local/...

You: Summarize what the pipeline does.
  [calling tool: summarize_nextflow_project(...)]
Agent: This is a 3-stage DSL2 pipeline: QC -> Assembly -> Annotation...

You: Is there a module missing for taxonomic classification?
  [calling tool: identify_missing_modules(...)]
Agent: Yes — the pipeline has QC and Assembly outputs but nothing consumes them for taxonomy. Here's a draft module design: ...

You: Use Kraken2 instead of the tool you picked, and write output to a TSV, not JSON.
  [calling tool: store_design_feedback(...)]
  [calling tool: update_module_design(...)]
Agent: Updated the design to use Kraken2 with TSV output. Want me to write it to disk?

You: Yes, write it to ./generated
  [calling tool: write_module_to_file(...)]
Agent: Wrote generated/modules/local/kraken2_classify/main.nf
```

Notes on this flow: the agent never writes files without you asking, it remembers your Kraken2/TSV preference in `feedback_memory.db` for future sessions (`recall_relevant_feedback`), and every read stayed inside `ALLOWED_ROOTS` / every write inside `MODULE_WRITE_ROOTS`.

---

## 7. Switching to a different model or provider

Because `chatOpenai.py` uses LangChain's `ChatOpenAI` client with a configurable `base_url`, any OpenAI-compatible chat-completions endpoint works — just change `BASE_URL`, `MODEL_NAME`, and `API_KEY` in `.env`:

```
# OpenAI directly
BASE_URL=https://api.openai.com/v1
MODEL_NAME=gpt-4o
API_KEY=sk-...

# PNNL AI Incubator gateway (default)
BASE_URL=https://ai-incubator-api.pnnl.gov
MODEL_NAME=gpt-5.4-project
API_KEY=sk-...
```

For a fully local/offline model instead (no external API calls at all), use `chatOllama.py` + `ollama.sh` (see `ollama.sh` for how to pull and serve a model with Apptainer/Ollama on an HPC node) — this is the path that *does* benefit from a SLURM/GPU allocation, since the model runs locally rather than on a remote gateway.

---

## 8. How the chat loop works (for the curious)

* Conversation history is kept in memory as a Python list of `SystemMessage` / `HumanMessage` / `AIMessage` / `ToolMessage` objects.
* Each turn sends the *full* history to the model so it has context.
* When the model wants to use a tool (e.g. `read_file`), LangChain returns a tool call instead of plain text; `chatOpenai.py` runs the corresponding Python function from `skills/`, feeds the result back to the model as a `ToolMessage`, and loops (up to 12 tool calls per turn) until the model gives a final plain-text answer.
* `/reset` throws away all history except the initial system prompt.

