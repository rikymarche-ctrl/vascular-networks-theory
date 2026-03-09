.PHONY: all paper1 paper2 clean setup help

PAPER1 = Beyond Murray's Law
PAPER2 = A Unified Variational Principle for Branching Transport Networks

help:
	@echo ""
	@echo "  make paper1   — compute + compile Paper 1"
	@echo "  make paper2   — compute + compile Paper 2"
	@echo "  make all      — both papers"
	@echo "  make setup    — install Python dependencies"
	@echo "  make clean    — remove LaTeX auxiliary files"
	@echo ""
	@echo "  Without make: python reproduce.py [paper1|paper2]"
	@echo ""

all: paper1 paper2

setup:
	pip install -r requirements.txt

paper1:
	python paper1-murray/compute_paper_murray.py
	cd paper1-murray/manuscript && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER1)" main.tex && \
	  bibtex "$(PAPER1)" && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER1)" main.tex && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER1)" main.tex

paper2:
	python paper2-variational/compute_paper_variational.py
	cd paper2-variational/manuscript && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER2)" main.tex && \
	  bibtex "$(PAPER2)" && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER2)" main.tex && \
	  pdflatex -interaction=nonstopmode -jobname "$(PAPER2)" main.tex

clean:
	find . \( -name "*.aux" -o -name "*.log" -o -name "*.bbl" -o -name "*.blg" \
	       -o -name "*.out" -o -name "*.toc" -o -name "*.fls" -o -name "*.fdb_latexmk" \
	       -o -name "*.synctex.gz" \) -delete
