check:
	Rscript -e 'devtools::check()'

test:
	Rscript -e 'devtools::test()'

manual: scCNAutils-manual.pdf

scCNAutils-manual.pdf: man/*.Rd
	rm scCNAutils-manual.pdf
	R CMD Rd2pdf . --no-preview -o scCNAutils-manual.pdf
