FROM continuumio/miniconda3:4.12.0

SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN conda create -n kalasanty_env -y \
    -c cheminfIBB -c anaconda -c openbabel \
    python=3.6 \
    h5py>=2.7 \
    openbabel=2.4.1 \
    tfbio=0.3 \
    scikit-image>=0.13 \
    numpy>=1.12 \
    scipy>=1 \
    keras=2.2.4 \
    tensorflow=1.11 \
    protobuf=3.6.1 \
    tqdm \
    && conda clean -afy

RUN git clone https://gitlab.com/cheminfIBB/kalasanty.git /opt/kalasanty \
    && conda run -n kalasanty_env pip install /opt/kalasanty

RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*

ENV KALASANTY_DIR=/opt/kalasanty
WORKDIR /workspace

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "kalasanty_env", "python"]
