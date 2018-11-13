# This file sets the location of configure.nems and modules.nems, and
# adds Make rules to create the tests/*.exe and tests/modules.* files.
# This file is included by the NEMS build system, within
# NEMS/GNUmakefile, just after platform logic is executed, but before
# the appbuilder file (if any) is read.

# IMPORTANT: This file MUST NOT CONTAIN any logic specific to building
# FV3, CCPP, FMS, or NEMS.  Otherwise, coupled FV3 applications will
# break.  There should only be logic specific to the NEMSfv3gfs test
# system and NEMSfv3gfs file naming in this makefile fragment.
#
# Logic specific to FV3, CCPP, FMS, or NEMS belong in NEMS/src/incmake.

# ----------------------------------------------------------------------
# Decide the conf and modulefile names.

CHOSEN_MODULE=$(BUILD_TARGET)/fv3

ifneq (,$(or $(findstring CCPP,$(COMPONENTS)),$(findstring CCPP,$(FV3_MAKEOPT)),$(findstring CCPP,$(FV3_MAKEOPT))))
  ifeq ($(CHOSEN_MODULE),theia.intel/fv3)
    override CHOSEN_MODULE=$(BUILD_TARGET)/fv3.intel-18.0.1.163
    $(warning Overriding CHOSEN_MODULE with $(CHOSEN_MODULE) because CCPP does not work on NOAA Theia with versions of the Intel compiler earlier than 18)
  endif
endif

CONFIGURE_NEMS_FILE=configure.fv3.$(BUILD_TARGET)

# ----------------------------------------------------------------------
# Copy the executable and modules.nems files into the tests/ directory
# if a TEST_BUILD_NAME is specified.

ifneq ($(TEST_BUILD_NAME),)
$(info Will copy modules.nems and NEMS.x as $(TEST_BUILD_NAME) under tests/)
$(ROOTDIR)/tests/$(TEST_BUILD_NAME).exe: $(NEMS_EXE)
	set -xe ; cp "$<" "$@"

$(ROOTDIR)/tests/modules.$(TEST_BUILD_NAME): $(NEMSDIR)/src/conf/modules.nems
	set -xe ; cp "$<" "$@"

configure: $(ROOTDIR)/tests/modules.$(TEST_BUILD_NAME) ;
build: $(ROOTDIR)/tests/$(TEST_BUILD_NAME).exe ;
endif

