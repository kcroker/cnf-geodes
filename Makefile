# Parameters used to generate the data!
alpha = 2.35
remnant_ml = 1
remnant_mh = 100
geode_k = 3

# File variables
data = data/0.05.dat data/0.2.dat data/1.dat data/2.dat
figures = figure-geode_loglog-crop.pdf ligo_population_figs

# Top level target is default
# If there is a bibliography and a document file, make the final
level3.pdf : $(figures) level3.tex tex_scratch/level3.bbl
	pdflatex  -interaction nonstopmode -output-directory=tex_scratch level3.tex
	pdflatex  -interaction nonstopmode -output-directory=tex_scratch level3.tex
	cp tex_scratch/level3.pdf ./

# ~~ FIGURES ~~
# Just do expicit targets for now
figure-geode_loglog-crop.pdf : $(data)
	gnuplot < figure-geode_loglog.gpt
	pdflatex  -interaction nonstopmode -output-directory=tex_scratch figure-geode_loglog.tex
	pdfcrop tex_scratch/figure-geode_loglog.pdf ./figure-geode_loglog-crop.pdf

# (Note that these will fail if you've not already prepared the data sets referenced in
# build_core_figures.sh)
ligo_population_figs: heatmap.gpt prototype.gpt build_core_figures.sh
	reprocess = $(file <ligo_figures_built)  
	$(if $(reprocess), echo "Skipping ligo figures", ./build_core_figures.sh && echo "built" > ligo_figures_built)

# ~~ DATA ~~
# Haha that actually worked
# the sed line extracts the redshift from the name of the data file (don't push your luck)
%.dat : $(fryer_distribution)
	./odb_analysis_revisited.py $(alpha) $(remnant_ml) $(remnant_mh) $(geode_k) $(shell echo $@ | sed 's/data\/\(.*\)\.dat/\1/') > $@

# If there is no aux, build once to get it
tex_scratch/level3.aux : level3.tex $(figures)
	pdflatex -interaction nonstopmode -output-directory=tex_scratch level3.tex

# ~~ HORSESHIT ~~
# If there is no bibliography, run bibtex on the aux
# bibtex is dumb and needs absolute paths?  Ayup.
tex_scratch/level3.bbl : tex_scratch/level3.aux
	BIBINPUTS="$(shell pwd)/tex_scratch/" TEXMFOUTPUT="tex_scratch/" bibtex tex_scratch/level3.aux

.PHONY : clean distclean init ligo_population_figs

# Initialize the working space
init :
	mkdir tex_scratch
	mkdir .gnuplot_scratch
	ln -s $(shell pwd)/thesis.bib tex_scratch/

# Clean up TeX products
# and link, because bibtex is dumb and can only search one directory...
clean :
	rm *.pdf
	rm tex_scratch/*
	rm .gnuplot_scratch/*
	ln -s $(shell pwd)/thesis.bib tex_scratch/

# In addition, remove the generated datasets!
dataclean : 
	rm $(data) $(fryer_distribution) data/*.tab || echo "End happily"
