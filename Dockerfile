FROM nmrlipids/gromacs:latest 

# Install Miniconda 
ARG MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bash \
      wget \
      bzip2 \
      git \
      libgomp1 \
      ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet ${MINICONDA_URL} -O /tmp/conda.sh && \
    bash /tmp/conda.sh -b -p /opt/conda && \
    rm /tmp/conda.sh && \
    /opt/conda/bin/conda clean --all --yes
ENV PATH="/opt/conda/bin:${PATH}"

# Copy requirements and install into Conda env
COPY Scripts/DatabankLib/requirements-dev.txt requirements-dev.txt
RUN conda create -n databank python=3.11 --yes && \
    conda run -n databank pip install --no-cache-dir -r requirements-dev.txt && \
    conda clean --all --yes && \
    rm requirements-dev.txt

# Expose Conda executable, Conda env, and GROMACS in PATH
ENV PATH="/opt/conda/bin:/opt/conda/envs/databank/bin:/usr/local/gromacs/bin:${PATH}"

# Create and switch to non-root runner user
RUN useradd -m -s /bin/bash runner
USER runner

# Enable Conda in interactive shells
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/runner/.bashrc && \
    echo "conda activate databank" >> /home/runner/.bashrc


# Set working directory and default command
WORKDIR /workspace
 