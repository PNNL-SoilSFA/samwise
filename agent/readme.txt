Files to run the SAMWISE AI Agent

BETA AI AGENT: 

If you would like to test out the AI Agent that can help you interrogate your genomes and their metabolisms, simply set up the OpenAI agent by running:

`bash setupAgentOpenai.sh`

This will generate a conda environment called `langchain-chat-openai` that you then need to activate with `conda activate langchain-chat-openai`. 

After this, you need to set up your environment file (see env.txt example in samwise-main/agent/) to include your API key as well as user settings. Once that is done, change the file name to .env instead of env.txt so that the agent can find it. 

Then, you can run the agent using `python chatOpenai.py` within the /agent/ folder (you need to cd /agent/ if you have not already).
