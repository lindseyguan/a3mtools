#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = a3mtools
PYTHON_VERSION = 3.10
PYTHON_INTERPRETER = python

#################################################################################
# COMMANDS                                                                      #
#################################################################################


## Install Python Dependencies
.PHONY: requirements
requirements:
	mamba env update --name $(PROJECT_NAME) --file environment.yml


## Delete all compiled Python files
.PHONY: clean
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8 and black (use `make format` to do formatting)
.PHONY: lint
lint:
	flake8 a3mtools
	isort --check --diff --profile black a3mtools
	black --check --config pyproject.toml a3mtools

## Format source code with black
.PHONY: format
format:
	black --config pyproject.toml a3mtools


## Set up python interpreter environment
.PHONY: create_environment
create_environment:
	mamba env create --name $(PROJECT_NAME) -f environment.yml
	
	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"
	

## Set up python interpreter environment
.PHONY: create_test_environment
create_test_environment:
	mamba env create --name "$(PROJECT_NAME)_test" -f ./devtools/conda-envs/test_env.yaml
	
	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)_test"


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################


## Make Dataset
.PHONY: data
data: requirements
	$(PYTHON_INTERPRETER) a3mtools/dataset.py


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys; \
lines = '\n'.join([line for line in sys.stdin]); \
matches = re.findall(r'\n## (.*)\n[\s\S]+?\n([a-zA-Z_-]+):', lines); \
print('Available rules:\n'); \
print('\n'.join(['{:25}{}'.format(*reversed(match)) for match in matches]))
endef
export PRINT_HELP_PYSCRIPT

help:
	@$(PYTHON_INTERPRETER) -c "${PRINT_HELP_PYSCRIPT}" < $(MAKEFILE_LIST)
