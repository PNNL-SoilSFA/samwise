"""
rag_skill.py
============

A simple "RAG skill" for LangChain: point it at a folder of mixed
documents (.txt, .md, .pdf, .docx) and it will:

  1. Load and chunk the files
  2. Embed the chunks and store them in a Chroma vector store on disk
  3. Expose a single LangChain `@tool` — `search_documents` — that the
     chat model can call to retrieve relevant chunks as context.

Re-running the script only rebuilds the index if the source folder has
changed since the last run (tracked via a simple manifest file), so
repeated runs are fast.

Configuration (via .env, in addition to the chat model's own vars):
    DOCS_DIR        - folder to index                (default: ./docs)
    VECTOR_DB_DIR   - where to persist the index      (default: ./chroma_db)
    EMBEDDING_MODEL - sentence-transformers model name
                       (default: all-MiniLM-L6-v2)
    RAG_TOP_K       - number of chunks to retrieve per query (default: 4)

Embeddings run locally via sentence-transformers, so this works
regardless of whether your chat-completions gateway also supports an
embeddings endpoint. If your gateway *does* support OpenAI-style
embeddings and you'd rather use those, swap `build_embeddings()` to
return a `langchain_openai.OpenAIEmbeddings` instance instead.
"""

import hashlib
import json
import os
from pathlib import Path

from langchain_chroma import Chroma
from langchain_community.document_loaders import (
    Docx2txtLoader,
    PyPDFLoader,
    TextLoader,
)
from langchain_core.tools import tool
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_text_splitters import RecursiveCharacterTextSplitter

DOCS_DIR = os.getenv("DOCS_DIR", "./docs")
VECTOR_DB_DIR = os.getenv("VECTOR_DB_DIR", "./chroma_db")
EMBEDDING_MODEL = os.getenv("EMBEDDING_MODEL", "all-MiniLM-L6-v2")
RAG_TOP_K = int(os.getenv("RAG_TOP_K", "4"))

MANIFEST_PATH = Path(VECTOR_DB_DIR) / "_manifest.json"

# Map file extensions to the loader that knows how to read them.
LOADERS_BY_EXT = {
    ".txt": TextLoader,
    ".md": TextLoader,
    ".pdf": PyPDFLoader,
    ".docx": Docx2txtLoader,
}


def build_embeddings() -> HuggingFaceEmbeddings:
    """Local embedding model — no extra API dependency required."""
    return HuggingFaceEmbeddings(model_name=EMBEDDING_MODEL)


def _scan_docs_dir(docs_dir: Path) -> dict:
    """
    Build a manifest of {relative_path: mtime} for every supported file
    in docs_dir, used to detect whether the index is stale.
    """
    manifest = {}
    for path in sorted(docs_dir.rglob("*")):
        if path.is_file() and path.suffix.lower() in LOADERS_BY_EXT:
            manifest[str(path.relative_to(docs_dir))] = path.stat().st_mtime
        else:
            manifest[str(path.relative_to(docs_dir))] = path.stat().st_mtime
    return manifest


def _load_documents(docs_dir: Path) -> list:
    """Load every supported file in docs_dir using the right loader."""
    documents = []
    for path in sorted(docs_dir.rglob("*")):
        if not path.is_file():
            continue
        loader_cls = LOADERS_BY_EXT.get(path.suffix.lower())
        if loader_cls is None:
            loader_cls = LOADERS_BY_EXT.get(".txt")
            continue
        try:
            loader = loader_cls(str(path))
            documents.extend(loader.load())
        except Exception as e:
            print(f"[rag_skill] Warning: failed to load {path}: {e}")
    return documents


def _index_is_stale(current_manifest: dict) -> bool:
    if not MANIFEST_PATH.exists():
        return True
    try:
        saved_manifest = json.loads(MANIFEST_PATH.read_text())
    except Exception:
        return True
    return saved_manifest != current_manifest


def build_or_load_vectorstore(force_rebuild: bool = False) -> Chroma:
    """
    Build a Chroma index from DOCS_DIR if it doesn't exist yet or the
    source files changed; otherwise load the existing persisted index.
    """
    docs_dir = Path(DOCS_DIR)
    docs_dir.mkdir(parents=True, exist_ok=True)
    Path(VECTOR_DB_DIR).mkdir(parents=True, exist_ok=True)

    embeddings = build_embeddings()
    current_manifest = _scan_docs_dir(docs_dir)

    needs_rebuild = force_rebuild or _index_is_stale(current_manifest)

    if not needs_rebuild:
        print(f"[rag_skill] Loading existing index from {VECTOR_DB_DIR}")
        return Chroma(
            persist_directory=VECTOR_DB_DIR,
            embedding_function=embeddings,
        )

    print(f"[rag_skill] (Re)building index from {docs_dir} ...")
    raw_docs = _load_documents(docs_dir)

    if not raw_docs:
        print(f"[rag_skill] Warning: no supported documents found in {docs_dir}. "
              f"Supported types: {', '.join(LOADERS_BY_EXT)}")

    splitter = RecursiveCharacterTextSplitter(
        chunk_size=1000,
        chunk_overlap=150,
    )
    chunks = splitter.split_documents(raw_docs)

    # Wipe any old persisted data, then create a fresh index.
    vectorstore = Chroma.from_documents(
        documents=chunks if chunks else [],
        embedding=embeddings,
        persist_directory=VECTOR_DB_DIR,
    )

    MANIFEST_PATH.write_text(json.dumps(current_manifest, indent=2))
    print(f"[rag_skill] Indexed {len(raw_docs)} document(s) -> {len(chunks)} chunk(s)")
    return vectorstore


# Module-level vectorstore, built once on import / first use.
_vectorstore: Chroma | None = None


def get_vectorstore() -> Chroma:
    global _vectorstore
    if _vectorstore is None:
        _vectorstore = build_or_load_vectorstore()
    return _vectorstore


@tool
def search_documents(query: str) -> str:
    """
    Search the indexed document collection for content relevant to the
    query and return the most relevant excerpts with their source files.
    Use this whenever the user asks a question that might be answered by
    the documents in the knowledge base (e.g. company docs, reports,
    notes) rather than general knowledge.
    """
    vectorstore = get_vectorstore()
    results = vectorstore.similarity_search(query, k=RAG_TOP_K)

    if not results:
        return "No relevant documents were found for this query."

    formatted = []
    for i, doc in enumerate(results, start=1):
        source = doc.metadata.get("source", "unknown source")
        formatted.append(f"[{i}] (source: {source})\n{doc.page_content}")

    return "\n\n".join(formatted)


if __name__ == "__main__":
    # Standalone usage: `python rag_skill.py` builds/refreshes the index
    # and runs a quick test query so you can sanity-check it on its own.
    get_vectorstore()
    test_query = "What is this document collection about?"
    print(f"\nTest query: {test_query}\n{'-' * 40}")
    print(search_documents.invoke(test_query))
