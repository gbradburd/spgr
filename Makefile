
%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux references.bib
	bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.eps : %.pdf
	inkscape $< --export-eps=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.pdf: %.asy
	asy -f pdf $<
