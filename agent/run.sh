#!/bin/bash

#conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
#conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
source ~/.bashrc
conda activate langchain-chat-openai

# Make sure the environment has every package the agent needs. This is
# a no-op (fast) if everything is already installed, so it's safe to
# leave in for every run -- it's what catches a stale/incomplete
# environment instead of failing later with "ModuleNotFoundError".
python -m pip install -r requirementsOpenai.txt

python chatOpenai.py

