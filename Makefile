template.pdf: template.tex main
	./main
	pdflatex template.tex
	pdflatex template.tex

template.tex: main

main: main.cpp
	cc main.cpp -o main

all: template.pdf

clean:
	rm template.* main

.PHONY: all clean test
