#/bin/bash

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda create -n langchain-chat-openai python=3.11 -y
source ~/.bashrc
conda activate langchain-chat-openai
pip install langchain langchain-openai python-dotenv
pip uninstall -y torch nvidia-nccl-cu13 nvidia-cublas nvidia-cudnn-cu13 nvidia-cusparselt-cu13 nvidia-nvshmem-cu13 nvidia-cuda-cupti nvidia-cuda-nvrtc nvidia-cuda-runtime nvidia-cufft nvidia-cufile nvidia-curand nvidia-cusolver nvidia-cusparse nvidia-nvtx nvidia-nvjitlink cuda-toolkit cuda-bindings cuda-pathfinder triton
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install --upgrade --force-reinstall --no-deps pillow

