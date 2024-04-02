FROM --platform=linux/amd64 mambaorg/micromamba:1.5.7
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"  # Not a preferred method!
USER root
RUN apt-get update && apt-get install -y procps && apt-get clean
USER $MAMBA_USER