CURDIR = $(shell pwd)

clean:
	rm -rf nest/nest_build;
	rm -rf nest/nest_module_build;
	cd nest/nest_src; make clean;

configure:
	cd nest/nest_src; \
	./bootstrap.sh; \
	./configure --prefix=$(CURDIR)/nest/nest_build

build:
	cd nest/nest_src; \
	make -j 4;

install:
	cd nest/nest_src; \
	make install;

submodule:
	cd nest/nest_module; ./bootstrap.sh;
	rm -rf nest/nest_module_build;
	mkdir nest/nest_module_build;
	cd nest/nest_module_build; \
	../nest_module/configure --with-nest=$(CURDIR)/nest/nest_build/bin/nest-config; \
	make install -j 4;

nest: configure build install submodule
	
python:
	pip install -r requirements.txt

all: nest python
