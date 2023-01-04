.PHONY: build

all: install path
install:
	pip3 install --user -r requirements.txt
path:
	@echo Installation completed.
	@echo
	@echo Now you can append the path to this directory to your PATH \
		variable to make the executable accesible from wherever.

