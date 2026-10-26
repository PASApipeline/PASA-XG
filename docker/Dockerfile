# PASA-XG image: minimap2, samtools, TransDecoder and PASA-XG on top of the LRAA image.
# LRAA_IMAGE must be an LRAA build that includes --pasa_mode and --preassembled_input
# (branch pasa_mode_v2 until it is released).
ARG LRAA_IMAGE=us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest
FROM ${LRAA_IMAGE}

ARG MAMBA_VERSION=1.5.8
RUN apt-get update && apt-get install -y --no-install-recommends curl bzip2 ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION} \
       | tar -xvj -C /usr/local bin/micromamba \
    && micromamba create -y -p /opt/pxg -c conda-forge -c bioconda \
       python=3.11 minimap2=2.28 samtools transdecoder=5.7.1 pysam pandas pyarrow biopython \
    && micromamba clean -a -y

COPY . /opt/PASA-XG
RUN /opt/pxg/bin/pip install --no-deps /opt/PASA-XG
ENV PATH=/opt/pxg/bin:$PATH
ENV PASA_XG_LRAA=LRAA
ENTRYPOINT ["pasa-xg"]
