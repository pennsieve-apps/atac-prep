.PHONY: help test fixture run clean validate

SERVICE_NAME ?= atac-prep
PYTHON ?= python3
export DOCKER_DEFAULT_PLATFORM = linux/amd64

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test      - pytest: the synthetic fixture through the processor (no Docker)"
	@echo "make fixture   - write the synthetic input and genome to data/fixture/"
	@echo "make run       - build the image and run it in Docker over data/fixture/linked/, writing to data/output/{atac,qc}/"
	@echo "make validate  - app.yml against the manifest schema (draft 2020-12), CITATION.cff with cffconvert (own venv)"
	@echo "make clean     - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture --layout both

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

validate:
	$(PYTHON) -c "import json, yaml, jsonschema; s = json.load(open('../workflow-single-cell/docs/app-manifest.v1.schema.json')); jsonschema.Draft202012Validator(s).validate(yaml.safe_load(open('app.yml'))); print('app.yml ok')"
	test -x .venv-cff/bin/cffconvert || (python3 -m venv .venv-cff && .venv-cff/bin/pip install -q cffconvert==2.0.0)
	.venv-cff/bin/cffconvert --validate && echo "CITATION.cff ok"

clean:
	rm -rf data/output
