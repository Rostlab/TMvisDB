ARG PYTHON_BASE=3.10-slim

FROM docker.io/python:$PYTHON_BASE AS builder

# Install PDM
RUN pip install -U pdm
# Disable update check
ENV PDM_CHECK_UPDATE=false

# Copy necessary files for dependency installation
COPY pyproject.toml pdm.lock README.md /project/
COPY src/ /project/src

# Set the working directory
WORKDIR /project

# Install dependencies and project into the local packages directory
RUN pdm install --check --prod --no-editable

FROM docker.io/python:$PYTHON_BASE
ARG TINI_VERSION="v0.19.0"
ARG STREAMLIT_PORT=8501
ARG GIT_HASH

COPY --from=builder /project/.venv/ /project/.venv
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini

RUN set -eux; \
	apt-get update; \
	apt-get install -y gosu; \
	rm -rf /var/lib/apt/lists/*; \
    chmod +x /tini;

ENV PATH="/project/.venv/bin:$PATH"
ENV GIT_HASH=${GIT_HASH:-dev}
ENV STREAMLIT_PORT=${STREAMLIT_PORT}
ENV DATABASE_URL="sqlite:///data/tmvis.db"
ENV MAINTENANCE_MODE="false"
ENV LOG_LEVEL="ERROR"

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

COPY src /project/src
COPY assets /project/assets
WORKDIR /project

EXPOSE ${STREAMLIT_PORT}

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["/tini", "--", "/entrypoint.sh"]
CMD ["sh", "-c", "streamlit run src/streamlitapp.py --browser.gatherUsageStats=false --server.port=${STREAMLIT_PORT} --server.address=0.0.0.0"]