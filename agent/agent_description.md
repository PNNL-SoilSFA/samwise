# SAMWISE Query Agent Description
Developer: Oceane Bel (obel@pnnl.gov)   
Project PIs: Josue Rodriguez-Ramos (josue.a.rodriguez.ramos@gmail.com), Kirsten Hofmockel (kirsten.hofmockel@pnnl.gov). 

## 1. Capabilities / limitations of the agent
 
__What it can genuinely do today__ (all tools currently wired into `ALL_TOOLS` in `chatOpenai.py`):
 
- __Filesystem browsing/reading__ (`fs_skill.py`): list/read files, find files by name pattern or by content, all constrained to `ALLOWED_ROOTS`.
- __Document RAG search__ (`rag_skill.py` + `fs_skill.py`'s `index_directory`): semantic search over `.txt/.md/.pdf/.docx` via a local Chroma vector store + HuggingFace embeddings.
- __Large file handling__ (`large_file_skill.py`): chunked/paginated reading and in-file search for files too big to load at once (with session state so it remembers where you left off).
- __Excel__ (`excell_skill.py`): streaming overview, search, and CSV/TSV export of `.xlsx/.xlsm` workbooks — designed to scale to very large sheets without blowing up memory.
- __R `.rds` files__ (`rds_skill.py`): reads via `pyreadr` (pure Python) or falls back to `Rscript` if R is installed — a very SAMWISE-relevant capability since tools like MicroTrait output `.rds`.
- __Python code understanding__ (`code_skill.py`): AST-based analysis, docstring/function extraction, project-level summaries.
- __Nextflow pipeline understanding__ (`nextflow_skill.py`): parses `.nf`/`nextflow.config` to explain processes, channels, parameters, profiles.
- __Nextflow module design__ (`module_design_skill.py` + `module_writer_skill.py`): detects gaps in a pipeline, drafts new module specs, iterates on your feedback, and — only on request — writes real `.nf` files, restricted to `MODULE_WRITE_ROOTS`.
- __Feedback memory__ (`feedback_memory_skill.py`): persists your corrections/preferences in a local SQLite DB and recalls them semantically in later sessions/turns.
 
__Real limitations, stated plainly:__
 
- __It cannot execute arbitrary analysis code right now.__ `exec_skill.py` (which supports running Python scripts/functions/snippets in a sandboxed subprocess with a timeout) exists in `skills/` but __is not imported or registered in `ALL_TOOLS`__ in `chatOpenai.py`. So today the agent can *read and explain* code/pipelines but cannot *run* them — I verified this directly against the tool list.
- __No general database or external API access.__ There's no SQL/Postgres/API-client skill; the only "databases" it touches are its own local SQLite state files (feedback memory, large-file reading state) and the Chroma vector index — both are local application state, not external data sources.
- __No live/experimental-metadata awareness beyond what's in files.__ It has no LIMS, ELN, or sample-tracking integration; any "metadata" comes only from whatever structured/text files you point it at (see below).
- __Everything is gated by allow-lists__ (`ALLOWED_ROOTS` for reads, `MODULE_WRITE_ROOTS` for writes) — by design, not a bug, but worth being explicit that it can't reach outside those paths at all.
- __No conversation persistence across process restarts__ — only feedback and RAG/file indexes persist; the raw chat transcript does not (no session save/resume file was found).
- __Single-turn tool budget cap of 12__ tool calls per user message (`max_tool_iterations=12` in `run_turn`) — very deep multi-step tasks could hit this and stop with a "too many tool-call iterations" message.
 
## 2. Can you give it a database similar to a SAMWISE output?
 
Yes, with caveats based on format:
 
- __Structured tabular output__ (`.rds` from R/MicroTrait, `.xlsx`, `.csv/.tsv/.json/.log`): fully supported via `rds_skill.py`, `excell_skill.py`, and `fs_skill.py`'s `read_file`/`find_files_by_content`. Just point `ALLOWED_ROOTS` at the output directory.
- __Prose-style reports/docs__ (`.txt/.md/.pdf/.docx`): supported via `rag_skill.py`/`index_directory` for semantic search.
- __Not supported directly__: actual database engines (SQLite/Postgres/MySQL dumps), HDF5, Parquet, BAM/VCF or other bioinformatics binary formats, or anything not covered by the loaders above — those would need a new skill.
 
## 3. Can it access external databases or call analysis tools?
 
__No, not currently.__ There is no skill for calling external databases (no SQL connector, no REST/API client skill) or invoking bioinformatics tools (BLAST, Kraken2, etc.) directly. As noted above, even generic script execution (`exec_skill.py`) is present in the codebase but not wired up — so right now the agent is strictly read/analyze/design, not "run this tool for me." Enabling `exec_skill.py` in `ALL_TOOLS` would be the fastest path to giving it the ability to *run* existing analysis scripts (still sandboxed to `ALLOWED_ROOTS`, with a timeout) — but that's a deliberate scope decision to make, not a bug fix, since it changes the trust model (it would then execute code, not just read it).
 
## 4. Can it see sample metadata (SAMWISE file structure / experimental conditions)?
 
Only indirectly, through the generic tools — there is __no SAMWISE-specific metadata schema__ baked into any skill. Concretely:
 
- If experimental conditions/sample metadata live in a Nextflow `params`/`nextflow.config`, `explain_nextflow_config` will surface them.
- If they live in a samplesheet/CSV/JSON/Excel file, `read_file`, `find_in_excel`, or `find_files_by_content` (grep-like content search) can retrieve them if you tell/ask the agent to look, or it discovers them while summarizing a pipeline.
- If they're written up in a report/README (prose), `search_documents`/RAG search will find them.
- There's no auto-discovery of "the metadata file" for a given SAMWISE run — the agent has to be pointed at, or stumble onto, the right file. Building a dedicated "read SAMWISE run metadata" skill (aware of a standard file/naming convention if one exists) is a natural, scoped future addition rather than something present now.