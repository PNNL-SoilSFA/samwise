Files to run the SAMWISE AI Agent

BETA AI AGENT: 

If you would like to test out the AI Agent that can help you interrogate your genomes and their metabolisms, simply set up the OpenAI agent by running:

`bash setupAgentOpenai.sh`

This will generate a conda environment called `langchain-chat-openai` that you then need to activate with `conda activate langchain-chat-openai`. 

After this, copy `.env.example` to `.env` in the `agent/` directory and set your API key and user settings there. The `.env` file is local configuration and must not be committed.

Then, you can run the agent using `python chatOpenai.py` within the /agent/ folder (you need to cd /agent/ if you have not already).
