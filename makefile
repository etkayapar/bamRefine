.PHONY: build

all: install build path
install:
	pip3 install --user -r requirements.txt
build:
	python3 setup.py build_ext --inplace
path:
	@echo Installation completed.
	@echo
	@echo Now you can append the path to this directory to your PATH \
		variable to make the executable accesible from wherever.

