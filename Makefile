PROJECT_DIR ?= $(CURDIR)
ENV_BIN := $(if $(BIOINFO_ENV_PREFIX),$(BIOINFO_ENV_PREFIX)/bin/,)
R_BIN ?= $(ENV_BIN)Rscript
PY_BIN ?= $(ENV_BIN)python

export R_BIN
export PY_BIN

.PHONY: all test qa reports environment public-qa manifest

all:
	bash scripts/run_all.sh "$(PROJECT_DIR)"

test:
	"$(R_BIN)" tests/run_tests.R --project "$(PROJECT_DIR)"

qa:
	"$(R_BIN)" R/05_scientific_qa.R --project "$(PROJECT_DIR)"

reports:
	"$(PY_BIN)" scripts/build_reports.py --project "$(PROJECT_DIR)"

environment:
	bash scripts/capture_environment.sh "$(PROJECT_DIR)"

public-qa:
	"$(PY_BIN)" scripts/validate_public_release.py --project "$(PROJECT_DIR)"

manifest:
	"$(PY_BIN)" scripts/finalize_release.py --project "$(PROJECT_DIR)"
