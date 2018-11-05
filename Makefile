check:
	Rscript -e 'devtools::check()'

test:
	Rscript -e 'devtools::test()'

manual: docs/scCNAutils-manual.pdf

docs/scCNAutils-manual.pdf: man/*.Rd
	rm -f docs/scCNAutils-manual.pdf
	R CMD Rd2pdf . --no-preview -o docs/scCNAutils-manual.pdf

flowcharts: docs/flowchart-cnasignal.svg

docs/flowchart-cnasignal.svg: docs/flowchart-cnasignal.mmd
	./node_modules/.bin/mmdc -t default -i docs/flowchart-cnasignal.mmd -o docs/flowchart-cnasignal.png

## install mermaid with 'npm install mermaid.cli'
