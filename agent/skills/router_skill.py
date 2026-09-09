# Add to skills/router_skill.py

"""
router_skill.py
===============

A helper tool that suggests which tools to use based on the user's question.
Call this first when unsure which tool is best for a query.
"""

import re
from langchain_core.tools import tool


@tool
def suggest_tools(query: str) -> str:
    """
    Analyzes a user question and suggests which tools to use.
    Call this when you're unsure which tool is best for a query.
    
    Args:
        query: The user's question or request
    
    Returns:
        Suggested tools and reasoning.
    """
    query_lower = query.lower()
    suggestions = []
    
    # Nextflow/Pipeline patterns
    pipeline_keywords = [
        "nextflow", "pipeline", "workflow", ".nf", "dram", "dram2",
        "nf-core", "process", "channel", "nextflow.config",
        "bioinformatics pipeline", "run the pipeline"
    ]
    if any(kw in query_lower for kw in pipeline_keywords):
        suggestions.append({
            "tools": ["find_nextflow_pipelines", "analyze_nextflow_pipeline", "summarize_nextflow_project"],
            "reason": "Query mentions Nextflow/pipeline concepts",
            "example": "summarize_nextflow_project('/path/to/pipeline')"
        })
    
    # Python code patterns
    python_keywords = [
        "python", ".py", "script", "function", "class", "import",
        "analyze code", "understand the code", "what does this script"
    ]
    if any(kw in query_lower for kw in python_keywords):
        suggestions.append({
            "tools": ["find_files (*.py)", "analyze_python_file", "explain_python_project"],
            "reason": "Query is about Python code",
            "example": "find_files('/path', '*.py') then analyze_python_file"
        })
    
    # R data patterns
    r_keywords = ["rds", ".rds", "r data", "readrds", "microtrait"]
    if any(kw in query_lower for kw in r_keywords):
        suggestions.append({
            "tools": ["read_rds_file"],
            "reason": "Query mentions R data files",
            "example": "read_rds_file('/path/to/file.rds')"
        })
    
    # File finding patterns
    find_keywords = [
        "find files", "look for", "search for", "where are",
        "list all", "show me the"
    ]
    if any(kw in query_lower for kw in find_keywords):
        suggestions.append({
            "tools": ["find_files", "find_files_by_content"],
            "reason": "Query is about finding/searching files",
            "example": "find_files('/path', '*.csv')"
        })
    
    # Document/knowledge patterns
    doc_keywords = [
        "documentation", "docs", "what do we know about",
        "search the knowledge", "in the documents"
    ]
    if any(kw in query_lower for kw in doc_keywords):
        suggestions.append({
            "tools": ["search_documents"],
            "reason": "Query is about indexed documentation",
            "example": "search_documents('query text')"
        })
    
    # Default suggestion
    if not suggestions:
        suggestions.append({
            "tools": ["find_files", "list_directory"],
            "reason": "No specific pattern detected - start with file discovery",
            "example": "find_files('/path', '*') or list_directory('/path')"
        })
    
    # Format output
    lines = ["Tool suggestions for this query:\n"]
    for i, sug in enumerate(suggestions, 1):
        lines.append(f"{i}. {', '.join(sug['tools'])}")
        lines.append(f"   Reason: {sug['reason']}")
        lines.append(f"   Example: {sug['example']}\n")
    
    return "\n".join(lines)
