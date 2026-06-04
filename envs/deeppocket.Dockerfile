FROM python:3.10-bullseye

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
        libgsl-dev \
        libxml2-dev \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/Discngine/fpocket.git /opt/fpocket \
    && sed -i 's|char \*residue_string\[|char residue_string[|g; s|strcpy(&residue_string,|strcpy(residue_string,|g' /opt/fpocket/src/fparams.c \
    && make -C /opt/fpocket \
    && make -C /opt/fpocket install \
    && rm -rf /opt/fpocket

RUN pip install --no-cache-dir \
        --index-url https://download.pytorch.org/whl/cpu \
        torch==2.3.0 \
        torchvision==0.18.0 \
    && pip install --no-cache-dir \
        biopython==1.79 \
        matplotlib==3.8.3 \
        matplotlib-inline==0.1.7 \
        molgrid==0.5.3 \
        numpy==1.23.5 \
        pluggy==1.4.0 \
        ProDy==2.4.1 \
        pyparsing==3.1.0 \
        scikit-image==0.23.2 \
        scikit-learn==1.5.0

RUN pip install --no-cache-dir wandb==0.17.0

RUN git clone https://github.com/devalab/DeepPocket.git /opt/DeepPocket \
    && python -c "from pathlib import Path; root=Path('/opt/DeepPocket'); p=root/'rank_pockets.py'; s=p.read_text(); s=s.replace('model.cuda()', 'model.to(\"cpu\")').replace(\"device='cuda'\", \"device='cpu'\").replace(\".to('cuda')\", \".to('cpu')\").replace('input_tensor[:,:14]', 'input_tensor[:, :14]'); p.write_text(s); p=root/'predict.py'; s=p.read_text(); s=s.replace('model.cuda()', 'model.to(\"cpu\")').replace('torch.load(args.class_checkpoint)', 'torch.load(args.class_checkpoint, map_location=torch.device(\"cpu\"))').replace('torch.load(args.seg_checkpoint)', 'torch.load(args.seg_checkpoint, map_location=torch.device(\"cpu\"))').replace('torch.cuda.empty_cache()', 'torch.cuda.empty_cache() if torch.cuda.is_available() else None'); p.write_text(s)"

ENV DEEPPOCKET_DIR=/opt/DeepPocket
WORKDIR /workspace

ENTRYPOINT ["python"]
