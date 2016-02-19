VERSION=2016.02
install: clean document
	cd ../ && R CMD INSTALL TnSeq

all: clean document vignette test reference check build install

reference:
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf

check:
	cd ../ && R CMD check TnSeq --no-build-vignettes

build:
	cd ../ && R CMD BUILD TnSeq

test:
	./run_tests.R

roxygen:
	rm -f NAMESPACE && Rscript -e "roxygen2::roxygenize()"

document:
	rm -f NAMESPACE && Rscript -e "devtools::document()"

vignette:
	Rscript -e "devtools::build_vignettes()"

clean_vignette:
	rm -f inst/doc/*

vt:	clean_vignette vignette install

clean:
	rm -rf TnSeq/
	rm -rf TnSeq.Rcheck/
	rm -rf TnSeq_${VERSION}.tar.gz

