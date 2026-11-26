# Build the manuscript. Figures are produced by code/paper_fig*.py into figs/paper/ and copied here unchanged.
FIGS = fig1_cuts fig2_octagons fig3_plus fig4_switch fig5_vacancy fig6_vertex
SRC  = fig1 fig2 fig3 fig4 fig5 fig6

all: figures main.pdf

figures:
	@mkdir -p figures
	@set -- $(SRC); for f in $(FIGS); do cp ../figs/paper/$$1.pdf figures/$$f.pdf; shift; done

main.pdf: main.tex refs.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

check: main.pdf
	@grep -c "undefined" main.log || true
	@grep -c "multiply defined" main.log || true
	@grep -c "Overfull" main.log || true
	@pdffonts main.pdf | awk 'NR>2 && ($$5!="yes" || $$2 ~ /Type 3/)' 

clean:
	latexmk -C main.tex

.PHONY: all figures check clean
