"""
feedback_memory_skill.py
========================

Persistent feedback memory for the chat agent. Stores user feedback in a
SQLite database with vector embeddings for semantic retrieval, allowing
the agent to learn from past interactions.

This skill enables:
  1. Storing feedback with context (what task/design it relates to)
  2. Retrieving relevant past feedback by semantic similarity
  3. Listing feedback history for specific topics or time ranges

Storage:
  - SQLite database for structured data (feedback text, timestamps, topics)
  - Vector embeddings stored alongside for semantic search
  - Persistent across sessions (survives restarts)

The combination of relational storage + vector search is a common pattern
known as "hybrid search" and is well-documented in the LangChain and
vector database literature.

References:
  - SQLite VSS extension: https://github.com/asg017/sqlite-vss
  - LangChain SQLite integration: https://python.langchain.com/docs/integrations/vectorstores/
  - Sentence Transformers: https://www.sbert.net/

Configuration (via .env):
    FEEDBACK_DB_PATH    - path to SQLite database (default: ./feedback_memory.db)
    FEEDBACK_EMBED_MODEL - embedding model name (default: all-MiniLM-L6-v2)
    FEEDBACK_TOP_K      - default number of results for similarity search (default: 5)

Dependencies:
    pip install sqlite-vss numpy
    (sentence-transformers already installed via rag_skill)
"""

import json
import os
import sqlite3
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
from langchain_core.tools import tool
from langchain_huggingface import HuggingFaceEmbeddings

# --- Configuration -----------------------------------------------------------

FEEDBACK_DB_PATH = Path(os.getenv("FEEDBACK_DB_PATH", "./feedback_memory.db"))
FEEDBACK_EMBED_MODEL = os.getenv("FEEDBACK_EMBED_MODEL", "all-MiniLM-L6-v2")
FEEDBACK_TOP_K = int(os.getenv("FEEDBACK_TOP_K", "5"))

# Embedding dimension for all-MiniLM-L6-v2 is 384
# Reference: https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2
EMBEDDING_DIM = 384


# --- Database Setup ----------------------------------------------------------

def _get_embeddings() -> HuggingFaceEmbeddings:
    """
    Build the embedding model. Uses the same approach as rag_skill.py
    for consistency.
    
    Reference: https://python.langchain.com/docs/integrations/text_embedding/huggingfacehub
    """
    return HuggingFaceEmbeddings(model_name=FEEDBACK_EMBED_MODEL)


# Module-level embedding model (loaded once)
_embeddings: HuggingFaceEmbeddings | None = None


def _get_embed_model() -> HuggingFaceEmbeddings:
    """Lazy-load the embedding model."""
    global _embeddings
    if _embeddings is None:
        _embeddings = _get_embeddings()
    return _embeddings


def _embed_text(text: str) -> np.ndarray:
    """
    Generate embedding vector for text.
    
    Returns a numpy array of shape (EMBEDDING_DIM,).
    """
    model = _get_embed_model()
    # HuggingFaceEmbeddings.embed_query returns a list of floats
    embedding = model.embed_query(text)
    return np.array(embedding, dtype=np.float32)


def _serialize_embedding(embedding: np.ndarray) -> bytes:
    """
    Serialize numpy array to bytes for SQLite BLOB storage.
    
    Using numpy's native binary format is efficient and well-documented:
    https://numpy.org/doc/stable/reference/generated/numpy.ndarray.tobytes.html
    """
    return embedding.tobytes()


def _deserialize_embedding(blob: bytes) -> np.ndarray:
    """Deserialize bytes back to numpy array."""
    return np.frombuffer(blob, dtype=np.float32)


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute cosine similarity between two vectors.
    
    Reference: Standard formula, see scikit-learn documentation
    https://scikit-learn.org/stable/modules/metrics.html#cosine-similarity
    """
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return float(dot_product / (norm_a * norm_b))


def _init_database() -> sqlite3.Connection:
    """
    Initialize the SQLite database with the feedback table.
    
    Schema design follows SQLite best practices:
    https://www.sqlite.org/lang_createtable.html
    
    The embedding is stored as a BLOB (binary) which is the standard
    approach for storing vectors in SQLite before using extensions.
    """
    FEEDBACK_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    
    conn = sqlite3.connect(str(FEEDBACK_DB_PATH))
    conn.row_factory = sqlite3.Row  # Enable column access by name
    
    # Create feedback table if it doesn't exist
    conn.execute("""
        CREATE TABLE IF NOT EXISTS feedback (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL,
            topic TEXT NOT NULL,
            subtopic TEXT,
            feedback_text TEXT NOT NULL,
            context TEXT,
            sentiment TEXT,
            source TEXT,
            embedding BLOB NOT NULL,
            metadata TEXT
        )
    """)
    
    # Create indexes for common queries
    # Reference: https://www.sqlite.org/lang_createindex.html
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_feedback_topic 
        ON feedback(topic)
    """)
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_feedback_timestamp 
        ON feedback(timestamp)
    """)
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_feedback_subtopic 
        ON feedback(subtopic)
    """)
    
    conn.commit()
    return conn


def _get_connection() -> sqlite3.Connection:
    """Get a database connection (creates DB if needed)."""
    return _init_database()


# --- Core Functions ----------------------------------------------------------

def _store_feedback(
    topic: str,
    feedback_text: str,
    subtopic: Optional[str] = None,
    context: Optional[str] = None,
    sentiment: Optional[str] = None,
    source: Optional[str] = None,
    metadata: Optional[dict] = None,
) -> int:
    """
    Store a feedback entry in the database.
    
    Returns the ID of the inserted row.
    """
    # Generate embedding from the combined text for better retrieval
    embed_text = f"{topic} {subtopic or ''} {feedback_text}"
    embedding = _embed_text(embed_text)
    
    conn = _get_connection()
    cursor = conn.execute(
        """
        INSERT INTO feedback 
        (timestamp, topic, subtopic, feedback_text, context, sentiment, source, embedding, metadata)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            datetime.now().isoformat(),
            topic,
            subtopic,
            feedback_text,
            context,
            sentiment,
            source,
            _serialize_embedding(embedding),
            json.dumps(metadata) if metadata else None,
        )
    )
    conn.commit()
    feedback_id = cursor.lastrowid
    conn.close()
    return feedback_id


def _search_feedback_by_similarity(
    query: str,
    topic_filter: Optional[str] = None,
    top_k: int = FEEDBACK_TOP_K,
) -> list[dict]:
    """
    Search feedback entries by semantic similarity.
    
    Uses brute-force cosine similarity search. For small to medium
    datasets (< 100k entries), this is efficient enough. For larger
    scale, consider sqlite-vss or a dedicated vector store.
    
    Reference on vector search approaches:
    https://www.pinecone.io/learn/vector-similarity/
    """
    query_embedding = _embed_text(query)
    
    conn = _get_connection()
    
    if topic_filter:
        rows = conn.execute(
            "SELECT * FROM feedback WHERE topic = ?",
            (topic_filter,)
        ).fetchall()
    else:
        rows = conn.execute("SELECT * FROM feedback").fetchall()
    
    conn.close()
    
    if not rows:
        return []
    
    # Compute similarities
    results = []
    for row in rows:
        stored_embedding = _deserialize_embedding(row["embedding"])
        similarity = _cosine_similarity(query_embedding, stored_embedding)
        results.append({
            "id": row["id"],
            "timestamp": row["timestamp"],
            "topic": row["topic"],
            "subtopic": row["subtopic"],
            "feedback_text": row["feedback_text"],
            "context": row["context"],
            "sentiment": row["sentiment"],
            "source": row["source"],
            "similarity": similarity,
            "metadata": json.loads(row["metadata"]) if row["metadata"] else None,
        })
    
    # Sort by similarity (descending) and return top_k
    results.sort(key=lambda x: x["similarity"], reverse=True)
    return results[:top_k]


def _get_feedback_by_topic(topic: str, subtopic: Optional[str] = None) -> list[dict]:
    """Retrieve all feedback for a specific topic (and optionally subtopic)."""
    conn = _get_connection()
    
    if subtopic:
        rows = conn.execute(
            "SELECT * FROM feedback WHERE topic = ? AND subtopic = ? ORDER BY timestamp DESC",
            (topic, subtopic)
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT * FROM feedback WHERE topic = ? ORDER BY timestamp DESC",
            (topic,)
        ).fetchall()
    
    conn.close()
    
    return [
        {
            "id": row["id"],
            "timestamp": row["timestamp"],
            "topic": row["topic"],
            "subtopic": row["subtopic"],
            "feedback_text": row["feedback_text"],
            "context": row["context"],
            "sentiment": row["sentiment"],
            "source": row["source"],
            "metadata": json.loads(row["metadata"]) if row["metadata"] else None,
        }
        for row in rows
    ]


def _list_topics() -> list[dict]:
    """List all unique topics with feedback counts."""
    conn = _get_connection()
    rows = conn.execute("""
        SELECT topic, COUNT(*) as count, MAX(timestamp) as last_updated
        FROM feedback
        GROUP BY topic
        ORDER BY count DESC
    """).fetchall()
    conn.close()
    
    return [
        {"topic": row["topic"], "count": row["count"], "last_updated": row["last_updated"]}
        for row in rows
    ]


# --- LangChain Tools ---------------------------------------------------------

@tool
def store_user_feedback(
    topic: str,
    feedback_text: str,
    subtopic: str = "",
    context: str = "",
    sentiment: str = "",
) -> str:
    """
    Store user feedback for future reference. Use this whenever the user
    provides feedback, corrections, preferences, or guidance that should
    be remembered for future interactions.
    
    Args:
        topic: Main category (e.g., "module_design", "code_style", 
               "pipeline_preference", "SORT_BAM")
        feedback_text: The actual feedback from the user
        subtopic: Optional sub-category (e.g., process name, file name)
        context: Optional context about what prompted this feedback
        sentiment: Optional: "positive", "negative", or "neutral"
    
    Returns:
        Confirmation message with the feedback ID.
    
    Examples:
        - User says "I prefer using biocontainers over quay.io"
          → store_user_feedback("container_preference", "Prefer biocontainers over quay.io")
        
        - User corrects a module design
          → store_user_feedback("module_design", "Output should emit both BAM and BAI", 
                                subtopic="SORT_BAM")
    """
    if not topic or not feedback_text:
        return "Error: Both 'topic' and 'feedback_text' are required."
    
    try:
        feedback_id = _store_feedback(
            topic=topic,
            feedback_text=feedback_text,
            subtopic=subtopic if subtopic else None,
            context=context if context else None,
            sentiment=sentiment if sentiment else None,
            source="user",
        )
        return (
            f"✓ Stored feedback #{feedback_id}\n"
            f"  Topic: {topic}" + (f" / {subtopic}" if subtopic else "") + "\n"
            f"  Feedback: {feedback_text[:100]}{'...' if len(feedback_text) > 100 else ''}"
        )
    except Exception as e:
        return f"Error storing feedback: {e}"


@tool
def recall_relevant_feedback(query: str, topic: str = "", top_k: int = 5) -> str:
    """
    Search for past user feedback relevant to the current task or question.
    Use this at the START of a task to check if the user has previously
    given relevant guidance, preferences, or corrections.
    
    Args:
        query: Description of what you're looking for (e.g., "container 
               preferences for bioinformatics tools", "feedback on SORT_BAM design")
        topic: Optional topic filter to narrow search
        top_k: Number of results to return (default: 5)
    
    Returns:
        Relevant past feedback entries with similarity scores.
    
    Examples:
        - Before designing a module:
          recall_relevant_feedback("module design preferences")
        
        - Before suggesting containers:
          recall_relevant_feedback("container preferences", topic="container_preference")
    """
    if not query:
        return "Error: 'query' is required."
    
    try:
        results = _search_feedback_by_similarity(
            query=query,
            topic_filter=topic if topic else None,
            top_k=top_k,
        )
        
        if not results:
            return f"No relevant feedback found for: '{query}'"
        
        lines = [f"Found {len(results)} relevant feedback entries:\n"]
        for i, r in enumerate(results, 1):
            sim_pct = f"{r['similarity'] * 100:.1f}%"
            lines.append(f"[{i}] (relevance: {sim_pct}) — {r['topic']}" + 
                        (f"/{r['subtopic']}" if r['subtopic'] else ""))
            lines.append(f"    {r['feedback_text']}")
            if r['context']:
                lines.append(f"    Context: {r['context'][:80]}...")
            lines.append(f"    (recorded: {r['timestamp'][:10]})")
            lines.append("")
        
        return "\n".join(lines)
    
    except Exception as e:
        return f"Error searching feedback: {e}"


@tool
def get_feedback_for_topic(topic: str, subtopic: str = "") -> str:
    """
    Retrieve ALL feedback entries for a specific topic. Use this when
    you need complete context about a particular subject, not just the
    most relevant entries.
    
    Args:
        topic: The topic to retrieve (e.g., "module_design", "SORT_BAM")
        subtopic: Optional subtopic filter
    
    Returns:
        All feedback entries for the topic, sorted by date (newest first).
    """
    if not topic:
        return "Error: 'topic' is required."
    
    try:
        results = _get_feedback_by_topic(topic, subtopic if subtopic else None)
        
        if not results:
            return f"No feedback found for topic: '{topic}'" + (f"/{subtopic}" if subtopic else "")
        
        lines = [f"Feedback for '{topic}'" + (f"/{subtopic}" if subtopic else "") + f" ({len(results)} entries):\n"]
        for r in results:
            lines.append(f"• [{r['timestamp'][:10]}] {r['feedback_text']}")
            if r['sentiment']:
                lines.append(f"  Sentiment: {r['sentiment']}")
            if r['context']:
                lines.append(f"  Context: {r['context'][:100]}...")
            lines.append("")
        
        return "\n".join(lines)
    
    except Exception as e:
        return f"Error retrieving feedback: {e}"


@tool
def list_feedback_topics() -> str:
    """
    List all topics that have stored feedback, with counts and last
    update times. Use this to see what feedback categories exist.
    
    Returns:
        List of topics with entry counts.
    """
    try:
        topics = _list_topics()
        
        if not topics:
            return "No feedback has been stored yet."
        
        lines = [f"Feedback topics ({len(topics)} total):\n"]
        for t in topics:
            lines.append(f"  • {t['topic']}: {t['count']} entries (last: {t['last_updated'][:10]})")
        
        return "\n".join(lines)
    
    except Exception as e:
        return f"Error listing topics: {e}"


@tool
def store_design_feedback(
    design_name: str,
    feedback_text: str,
    feedback_type: str = "general",
    approved: bool = False,
) -> str:
    """
    Store feedback specifically about a module design. This is a
    convenience wrapper around store_user_feedback that's optimized
    for the module design workflow.
    
    Args:
        design_name: Name of the design (e.g., "SORT_BAM", "INDEX_BAM")
        feedback_text: The user's feedback
        feedback_type: Type of feedback: "inputs", "outputs", "resources", 
                      "container", "script", "general", "approval"
        approved: Whether this feedback indicates approval
    
    Returns:
        Confirmation message.
    """
    if not design_name or not feedback_text:
        return "Error: Both 'design_name' and 'feedback_text' are required."
    
    sentiment = "positive" if approved else None
    if approved:
        feedback_type = "approval"
    
    try:
        feedback_id = _store_feedback(
            topic="module_design",
            subtopic=design_name,
            feedback_text=feedback_text,
            context=f"Feedback type: {feedback_type}",
            sentiment=sentiment,
            source="user",
            metadata={"feedback_type": feedback_type, "approved": approved},
        )
        
        status = "✅ APPROVED" if approved else "📝 Feedback recorded"
        return (
            f"{status} for design '{design_name}' (#{feedback_id})\n"
            f"  Type: {feedback_type}\n"
            f"  Feedback: {feedback_text[:150]}{'...' if len(feedback_text) > 150 else ''}"
        )
    except Exception as e:
        return f"Error storing design feedback: {e}"


# --- Optional: Database maintenance tools ------------------------------------

@tool
def get_feedback_stats() -> str:
    """
    Get statistics about the feedback database. Useful for understanding
    how much feedback has been collected.
    
    Returns:
        Database statistics including total entries, topics, date range.
    """
    try:
        conn = _get_connection()
        
        total = conn.execute("SELECT COUNT(*) FROM feedback").fetchone()[0]
        topics = conn.execute("SELECT COUNT(DISTINCT topic) FROM feedback").fetchone()[0]
        
        if total > 0:
            oldest = conn.execute("SELECT MIN(timestamp) FROM feedback").fetchone()[0]
            newest = conn.execute("SELECT MAX(timestamp) FROM feedback").fetchone()[0]
            
            sentiments = conn.execute("""
                SELECT sentiment, COUNT(*) 
                FROM feedback 
                WHERE sentiment IS NOT NULL 
                GROUP BY sentiment
            """).fetchall()
        else:
            oldest = newest = None
            sentiments = []
        
        conn.close()
        
        lines = ["Feedback Database Statistics:\n"]
        lines.append(f"  Total entries: {total}")
        lines.append(f"  Unique topics: {topics}")
        
        if total > 0:
            lines.append(f"  Date range: {oldest[:10]} to {newest[:10]}")
            if sentiments:
                lines.append(f"  Sentiments: " + ", ".join(f"{s[0]}={s[1]}" for s in sentiments))
        
        lines.append(f"\n  Database location: {FEEDBACK_DB_PATH}")
        
        return "\n".join(lines)
    
    except Exception as e:
        return f"Error getting stats: {e}"
