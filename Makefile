# Parameters used to generate the data!
alpha = 2.35
remnant_ml = 1
remnant_mh = 100
geode_k = 3

# File variables
data = data/0.05.dat data/0.2.dat data/1.dat data/2.dat
figures = figure-geode_loglog-crop.pdf

# Top level target is default

# ~~ FIGURES ~~
# Just do expicit targets for now
figure-geode_loglog-crop.pdf : $(data)
	gnuplot < figure-geode_loglog.gpt
	pdflatex  -interaction nonstopmode -output-directory=tex_scratch figure-geode_loglog.tex
	pdfcrop tex_scratch/figure-geode_loglog.pdf ./figure-geode_loglog-crop.pdf

# ~~ DATA ~~
# Haha that actually worked
# the sed line extracts the redshift from the name of the data file (don't push your luck)
%.dat :
	./odb_analysis_revisited.py $(alpha) $(remnant_ml) $(remnant_mh) $(geode_k) $(shell echo $@ | sed 's/data\/\(.*\)\.dat/\1/') > $@

.PHONY : clean distclean init 

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
