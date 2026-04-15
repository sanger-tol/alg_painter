# ---- Build stage ----
FROM python:3.12-slim AS builder

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# Install dependencies first (cached layer)
COPY pyproject.toml uv.lock .python-version ./
RUN uv sync --frozen --no-install-project --no-dev

# Copy source and install the project itself
COPY src/ src/
COPY README.md ./
RUN uv sync --frozen --no-dev --no-editable

# ---- Runtime stage ----
FROM python:3.12-slim AS runtime

WORKDIR /app

# Copy the virtual environment from builder
COPY --from=builder /app/.venv /app/.venv

# Ensure the venv is on PATH
ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONUNBUFFERED=1

# Matplotlib needs a writable config directory
ENV MPLCONFIGDIR=/tmp/matplotlib
