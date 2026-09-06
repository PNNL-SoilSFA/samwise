"""
Example: An interactive terminal chat window using LangChain with an
OpenAI-compatible model, parameterized via environment variables loaded
from a .env file.

Expects a .env file in the same directory with:
    BASE_URL=https://ai-incubator-api.pnnl.gov
    MODEL_NAME=gpt-5.4-project
    API_KEY=sk-*

Usage:
    python main.py

Chat commands:
    /reset   - clear conversation history and start over
    /exit    - quit the chat (also: /quit, Ctrl+C, Ctrl+D)
"""

import os
import sys

from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage, ToolMessage

from skills.rag_skill import search_documents
from skills.fs_skill import list_directory, read_file, index_directory, search_indexed_directory
from skills.rds_skill import read_rds_file
from skills.code_skill import analyze_python_file, explain_python_project, get_function_source 
from skills.fs_skill import find_files, find_files_by_content, list_directory, read_file, index_directory, search_indexed_directory

from skills.nextflow_skill import (
    find_nextflow_pipelines,
    analyze_nextflow_pipeline,
    explain_nextflow_config,
    summarize_nextflow_project,
    get_nextflow_process,
    explain_pipeline_directory,  # NEW - one-shot explainer
)

from skills.feedback_memory_skill import (
    store_user_feedback,
    recall_relevant_feedback,
    get_feedback_for_topic,
    list_feedback_topics,
    store_design_feedback,
    get_feedback_stats,
)

# Nextflow module design: identify missing modules, design them, and
# refine the designs using user feedback (before any implementation).
from skills.module_design_skill import (
    identify_missing_modules,
    gather_module_context,
    save_module_design,
    update_module_design,
    get_module_design,
    list_module_designs,
    render_module_design,
    record_design_feedback,
)

from skills.large_file_skill import (
    open_large_file,
    read_next_chunk,
    read_file_chunk,
    get_file_reading_status,
    reset_file_reading,
    read_file_range,
    search_in_large_file,
    get_file_summary,
    list_reading_sessions,
)

from skills.excell_skill import (
    read_excel, 
    find_in_excel, 
    export_excel_sheet,
)

from skills.module_writer_skill import (
    write_module_to_file,
    write_draft_module,
    write_all_approved_modules,
    preview_module_file,
    create_pipeline_scaffold,
    list_written_modules,
    show_write_permissions,
)

from skills.router_skill import suggest_tools

ALL_TOOLS = [
    # Router
    suggest_tools,

    # RAG / document search
    search_documents,

    # Filesystem
    find_files,
    find_files_by_content,
    list_directory,
    read_file,
    index_directory,
    search_indexed_directory,

    # R data files
    read_rds_file,

    # Python code analysis
    analyze_python_file,
    explain_python_project,
    get_function_source,

    # Nextflow pipelines
    explain_pipeline_directory,  # 🎯 One-shot, avoids iteration
    find_nextflow_pipelines,
    analyze_nextflow_pipeline,
    explain_nextflow_config,
    summarize_nextflow_project,
    get_nextflow_process,

    # Nextflow module design (identify missing → design → get feedback)
    identify_missing_modules,
    gather_module_context,
    save_module_design,
    update_module_design,
    get_module_design,
    list_module_designs,
    render_module_design,
    record_design_feedback,

    # Feedback memory
    store_user_feedback,
    recall_relevant_feedback,
    get_feedback_for_topic,
    list_feedback_topics,
    store_design_feedback,
    get_feedback_stats,

    # Large file reading
    open_large_file,
    read_next_chunk,
    read_file_chunk,
    get_file_reading_status,
    reset_file_reading,
    read_file_range,
    search_in_large_file,
    get_file_summary,
    list_reading_sessions,

    read_excel,
    find_in_excel,
    export_excel_sheet,

    write_module_to_file,
    write_draft_module,
    write_all_approved_modules,
    preview_module_file,
    create_pipeline_scaffold,
    list_written_modules,
    show_write_permissions,
]

TOOLS_BY_NAME = {tool.name: tool for tool in ALL_TOOLS}
# --- Load environment variables from .env -----------------------------------
load_dotenv()  # looks for a .env file in the current working directory

BASE_URL = os.getenv("BASE_URL")
MODEL_NAME = os.getenv("MODEL_NAME")
API_KEY = os.getenv("API_KEY")

missing = [name for name, val in
           [("BASE_URL", BASE_URL), ("MODEL_NAME", MODEL_NAME), ("API_KEY", API_KEY)]
           if not val]
if missing:
    sys.exit(f"Missing required environment variable(s): {', '.join(missing)}. "
              f"Check your .env file.")


def extract_text(content) -> str:
    """
    Normalize a message/chunk's `.content` into plain text.

    Most OpenAI-compatible backends return `content` as a plain string,
    but some return a list of content-block dicts, e.g.:
        [{'type': 'text', 'text': 'Hello', ...}, ...]
    This handles both shapes (and skips non-text blocks/empty fragments)
    so callers always get a plain string back.
    """
    if content is None:
        return ""
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts = []
        for block in content:
            if isinstance(block, str):
                parts.append(block)
            elif isinstance(block, dict) and block.get("type") == "text":
                parts.append(block.get("text", ""))
        return "".join(parts)
    return str(content)


def build_chat_model() -> ChatOpenAI:
    """
    Build a ChatOpenAI client pointed at a custom OpenAI-compatible endpoint.

    `base_url` is what makes this work against something other than
    api.openai.com — e.g. an internal gateway/incubator API that speaks
    the OpenAI chat-completions protocol.
    """
    llm = ChatOpenAI(
        model=MODEL_NAME,
        api_key=API_KEY,
        base_url=BASE_URL,
        temperature=0.7,
        # max_tokens=1024,   # uncomment / adjust as needed
        # timeout=60,        # uncomment / adjust as needed
    )
    # Give the model access to the RAG retrieval tool. It decides on its
    # own, per message, whether a question needs document context.
    #return llm.bind_tools([search_documents])
    return llm.bind_tools(ALL_TOOLS)


#TOOLS_BY_NAME = {tool.name: tool for tool in [search_documents]}


def run_turn(llm: ChatOpenAI, history: list, max_tool_iterations: int = 12) -> AIMessage:
    """
    Send `history` to the model and handle any tool calls it makes.

    The model may respond with one or more tool calls instead of a final
    answer (e.g. it decides it needs document context). We execute each
    requested tool, append the results as ToolMessages, and ask the model
    again — repeating until it returns a plain text answer or we hit
    `max_tool_iterations` as a safety cap against infinite loops.

    Note: the module-design workflow can chain several tool calls in a
    single turn (identify → gather context per candidate → save), so the
    cap is set generously. The model is instructed to end its turn (a
    plain text answer) once it has a design to present, so it can collect
    user feedback before doing anything further.
    """
    for _ in range(max_tool_iterations):
        ai_message = llm.invoke(history)
        history.append(ai_message)

        if not ai_message.tool_calls:
            return ai_message  # final answer, no further tool calls needed

        for call in ai_message.tool_calls:
            tool_fn = TOOLS_BY_NAME.get(call["name"])
            if tool_fn is None:
                result = f"Error: unknown tool '{call['name']}'"
            else:
                print(f"  [calling tool: {call['name']}({call['args']})]")
                try:
                    result = tool_fn.invoke(call["args"])
                except Exception as e:
                    result = f"Error running tool '{call['name']}': {e}"

            history.append(ToolMessage(content=str(result), tool_call_id=call["id"]))

    # Safety net: too many tool-call rounds without a final answer.
    fallback = AIMessage(content="(Stopped after too many tool-call iterations.)")
    history.append(fallback)
    return fallback


def main() -> None:
    llm = build_chat_model()

    # Conversation history. The SystemMessage is kept at index 0 and
    # persists across the session unless the user resets it.
    #system_prompt = (
    #    "You are a concise, helpful assistant. You have access to a "
    #    "'search_documents' tool that retrieves relevant excerpts from an "
    #    "internal document collection. Use it when the user's question "
    #    "could be answered by those documents; otherwise answer directly."
    #)
    #system_prompt = (
    #    "You are a concise, helpful assistant. You have these tools:\n"
    #    "- search_documents: pre-indexed docs from DOCS_DIR\n"
    #    "- list_directory / read_file: browse and read files in allow-listed "
    #    "directories (e.g. project output folders) the user names\n"
    #    "- index_directory / search_indexed_directory: build/query a semantic "
    #    "index over prose documents (.txt/.md/.pdf/.docx) in an arbitrary "
    #    "allow-listed directory\n"
    #    "For structured output like CSV/TSV/JSON, prefer list_directory + "
    #    "read_file over indexing."
    #)
    #system_prompt = """You are a helpful research assistant with filesystem, code, and pipeline analysis tools.

    system_prompt = """You are a helpful research assistant with filesystem, code, and pipeline analysis tools.

IMPORTANT: FEEDBACK MEMORY
At the START of any design or significant task, use recall_relevant_feedback() to check
if the user has previously given relevant guidance. Apply past feedback to current work.

When the user provides feedback, corrections, or preferences:
→ store_user_feedback(topic, feedback_text, subtopic, context)
→ For module designs specifically: store_design_feedback(design_name, feedback_text, feedback_type)

TOOL SELECTION (follow this order):

1. **Nextflow/Pipelines** (DRAM, DRAM2, nf-core, any .nf files):
   → summarize_nextflow_project(path) for overview
   → analyze_nextflow_pipeline(path) for specific .nf file
   → find_nextflow_pipelines(path) to discover pipelines

2. **Python code** (.py files, scripts):
   → find_files(path, "*.py") to find scripts
   → analyze_python_file(path) to understand code

3. **R data** (.rds files):
   → read_rds_file(path)

4. **Finding files**:
   → find_files(path, pattern) - PREFERRED, fast
   → list_directory(path) - only for full listings

5. **Documentation search** (only for indexed docs/reports):
   → search_documents(query)

6. **Designing missing Nextflow modules** (when the user asks to find,
   design, or add missing modules/processes to a pipeline):
   → identify_missing_modules(path)      # FIRST — read-only static scan
   → gather_module_context(path, NAME)   # for each candidate module
   → draft a spec from that context, then save_module_design(design)

   Then STOP and present each design to the user. Ask explicitly whether
   the inputs/outputs, resource directives, and container/conda match
   their intent. Do NOT call any more tools until they respond — end your
   turn with a plain-text answer so you can collect their feedback.

   When the user gives feedback:
   → get_module_design(NAME) to reload the design if needed
   → update_module_design(NAME, changes, feedback) to apply their changes
     (pass their comment in `feedback`); then show the revised draft and
     confirm it now matches their intent.

   Use record_design_feedback(NAME, ..., approved=True) to mark a design
   approved ONLY once the user is happy with it. NEVER write .nf files
   into the pipeline — implementation is a separate step that only
   happens after a design is approved.

7. **Reading large files** (logs, large source files, data files):
   → get_file_summary(path) to see file size and preview
   → open_large_file(path) then read_next_chunk(path) repeatedly for sequential reading
   → read_file_chunk(path, chunk_number=N) for random access
   → search_in_large_file(path, term) to find specific content
   → read_file_range(path, start_line, end_line) for specific line ranges

   For small files (< 8KB), continue using read_file() from fs_skill.

8. **Writing designed modules to disk** (when user asks to save/write modules):
   → preview_module_file(name) to review before writing
   → write_module_to_file(name, output_dir) for approved modules
   → write_draft_module(name, output_dir) for drafts (with warning)
   → write_all_approved_modules(output_dir) to write all at once
   → create_pipeline_scaffold(output_dir, name) to create a new pipeline structure

   IMPORTANT: Only write approved designs by default. If user wants to write
   a draft, use write_draft_module() which adds a warning header.
   
   Write permissions are restricted to MODULE_WRITE_ROOTS directories.
   Use show_write_permissions() to see allowed locations.

When user asks to "explain" a pipeline like DRAM2:
→ First use find_nextflow_pipelines to locate it
→ Then use summarize_nextflow_project or analyze_nextflow_pipeline

Do NOT use search_documents for code or pipeline questions.
"""
    history = [SystemMessage(content=system_prompt)]

    print("=" * 60)
    print(f" Chat with {MODEL_NAME}  (via {BASE_URL})")
    print(" Type your message and press Enter.")
    print(" Commands: /reset to clear history, /exit to quit.")
    print(f" RAG tool active — documents indexed from: {os.getenv('DOCS_DIR', './docs')}")
    print("=" * 60)

    while True:
        try:
            user_input = input("\nYou: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye!")
            break

        if not user_input:
            continue

        if user_input.lower() in ("/exit", "/quit"):
            print("Goodbye!")
            break

        if user_input.lower() == "/reset":
            history = [SystemMessage(content=system_prompt)]
            print("(conversation history cleared)")
            continue

        history.append(HumanMessage(content=user_input))
        snapshot_len = len(history)  # for rollback if this turn fails

        try:
            response = run_turn(llm, history)
        except Exception as e:
            print(f"\n[Error contacting model: {e}]")
            del history[snapshot_len - 1:]  # drop this turn entirely (user msg + any partial tool exchange)
            continue

        print(f"Assistant: {extract_text(response.content)}")


if __name__ == "__main__":
    main()
