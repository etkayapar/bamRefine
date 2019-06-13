neogene_install:
	sbatch install_dependencies.slurm

baggins_install:
	bash install_dependencies.sh

baggins_build:
	python3 setup.py build_ext --inplace

