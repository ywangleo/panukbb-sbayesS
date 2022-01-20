FROM ubuntu:latest
RUN apt-get update && apt-get install -y \
    python3.8 \
    wget \
    unzip \
    && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir sbayesS && \ 
    cd sbayesS && \
    wget https://cnsgenomics.com/software/gctb/download/gctb_2.03beta_Linux.zip && \
    unzip gctb_2.03beta_Linux.zip && \
    rm -rf gctb_2.03beta_Linux.zip

WORKDIR /sbayesS

COPY format_to_ma.py /sbayesS
