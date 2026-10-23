PYTHON ?= python3

.PHONY: setup setup-r setup-python validate main main-with-tree supplementary figure5 manifest check-manifest

setup: setup-r setup-python

setup-r:
	Rscript code/setup/install_r_dependencies.R

setup-python:
	$(PYTHON) -m pip install -r requirements.txt

validate:
	Rscript -e 'fs <- list.files("code", pattern="[.]R$$", recursive=TRUE, full.names=TRUE); for (f in fs) parse(file=f); cat("Parsed", length(fs), "R scripts\n")'
	$(PYTHON) -m compileall -q code
	Rscript code/validate_binary_inputs.R
	$(PYTHON) code/validate_repository.py

main:
	Rscript run_all.R

main-with-tree:
	Rscript run_all.R --include-tree

supplementary:
	Rscript run_supplementary.R

figure5:
	$(PYTHON) code/assemble_figure5_with_family.py

manifest:
	find . -path './.git' -prune -o -path '*/__pycache__' -prune -o -type f ! -name 'MANIFEST_SHA256.txt' ! -name '.DS_Store' ! -name '*.pyc' -print0 | LC_ALL=C sort -z | xargs -0 shasum -a 256 > MANIFEST_SHA256.txt

check-manifest:
	shasum -a 256 -c MANIFEST_SHA256.txt
