.DEFAULT_GOAL := help

PYTHON ?= python3
PAGES_BIND ?= 127.0.0.1
PAGES_PORT ?= 8000

.PHONY: help workflow build host clean

help:
	@printf '%s\n' \
	  'make build     Build the complete GitHub Pages site locally.' \
	  'make host     Serve the built site at http://$(PAGES_BIND):$(PAGES_PORT).' \
	  'make clean     Remove the generated site and local build workspace.' \
	  'Override NCAR_PHYSICS_REF, V7_PHYSICS_REF, or UFS_PHYSICS_REF to select refs.' \
	  'See README.md for prerequisites and configuration.'

NCAR_PHYSICS_REF ?= main
V7_PHYSICS_REF ?= v7.0.1
UFS_PHYSICS_REF ?= release/srw-v3
build:
	"$(PYTHON)" build_doxygen.py \
	  --ncar-physics-ref "$(NCAR_PHYSICS_REF)" \
	  --v7-physics-ref "$(V7_PHYSICS_REF)" \
	  --ufs-physics-ref "$(UFS_PHYSICS_REF)"

host:
	@test -f .build/doxygen/site/index.html || { \
	  printf '%s\n' 'No built site found. Run make build first.' >&2; \
	  exit 1; \
	}
	@printf '%s\n' 'Serving http://$(PAGES_BIND):$(PAGES_PORT) (Ctrl-C to stop).'
	"$(PYTHON)" -m http.server "$(PAGES_PORT)" --bind "$(PAGES_BIND)" --directory .build/doxygen/site

clean:
	rm -rf .build/doxygen
	@printf '%s\n' 'Removed tools/doxygen/.build/doxygen.'
