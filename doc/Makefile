LFLAGS  = -interaction=errorstopmode -halt-on-error
PROJECT	= hspf_ros_documentation
FILEDIR = files/
SRC     = $(PROJECT).tex

FILES = $(PROJECT).aux \
	$(PROJECT).bbl \
	$(PROJECT).blg \
	$(PROJECT).loa \
	$(PROJECT).lof \
	$(PROJECT).log \
	$(PROJECT).lol \
	$(PROJECT).lot \
	$(PROJECT).nav \
	$(PROJECT).out \
	$(PROJECT).snm \
	$(PROJECT).toc \

default:
	pdflatex $(LFLAGS) $(SRC)
	bibtex $(PROJECT)
	pdflatex $(LFLAGS) $(SRC)
	pdflatex $(LFLAGS) $(SRC)
	evince $(PROJECT).pdf

.PHONY: auto

auto:
	pdflatex $(LFLAGS) $(SRC)

.PHONY: clean

clean:
	@for file in $(FILES); \
		do \
		rm -f $$file; \
		done
