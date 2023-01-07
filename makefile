.PHONY: build

all:
devenv:
	pip3 install --user -r requirements-dev.txt
build: devenv
	python3 -m build
install:
	pip3 install --user -e .
test: devenv install
	python3 -m pytest
publish-test: test build
	python3 -m twine upload --repository testpypi dist/*
publish: test build
	python3 -m twine upload dist/*
clean:
	rm -rf dist/*
	pip3 uninstall -y bamrefine
