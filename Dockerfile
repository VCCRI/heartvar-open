ARG VEP_RELEASE=release_113.0
FROM ensemblorg/ensembl-vep:${VEP_RELEASE}

USER root
WORKDIR /app

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    PATH="/opt/venv/bin:/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/htslib:${PATH}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        software-properties-common gnupg ca-certificates \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3.11-dev \
        gcc libcurl4-openssl-dev zlib1g-dev \
    && python3.11 -m venv /opt/venv \
    && /opt/venv/bin/pip install --upgrade pip \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN /opt/venv/bin/pip install -r requirements.txt

RUN /opt/venv/bin/pip install --no-deps hgvs==1.5.7 cdot==0.2.31 \
    && /opt/venv/bin/python -c "import hgvs.parser, cdot.hgvs.dataproviders, psycopg2; print('hgvs stack OK')"

COPY backend/ ./backend/
COPY static/ ./static/
COPY index.html ./

RUN useradd --create-home --uid 10001 heartvar \
    && chown -R heartvar:heartvar /app
USER heartvar

EXPOSE 8000

ENTRYPOINT []

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD /opt/venv/bin/python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8000/health').status==200 else 1)"

CMD ["/opt/venv/bin/uvicorn", "backend.app:app", "--host", "0.0.0.0", "--port", "8000"]
