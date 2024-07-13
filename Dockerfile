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

COPY src /project/src

RUN useradd -m user && chown -R user:user /project


COPY --from=builder /project/.venv/ /project/.venv
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

ENV PATH="/project/.venv/bin:$PATH"
ENV GIT_HASH=${GIT_HASH:-dev}
ENV STREAMLIT_PORT=${STREAMLIT_PORT}
ENV DATABASE_URL="sqlite:///data/tmvis.db"
ENV MAINTENANCE_MODE="false"

WORKDIR /project

EXPOSE ${STREAMLIT_PORT}

USER user
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health
CMD ["/tini", "--", "sh", "-c", "streamlit run src/streamlitapp.py --server.port=${STREAMLIT_PORT} --server.address=0.0.0.0"]