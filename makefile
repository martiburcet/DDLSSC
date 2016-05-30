#=======================================================================
# Rules to build libraries and programs
.PHONY: default all debug release 
#config
.PHONY: clean clean_all clean_release clean_debug

default: debug

all: debug release

#config:
#	@make --no-print-directory -f make.drivers config

#newfile:
#	@make --no-print-directory -f make.drivers newfile

debug:
	@make --no-print-directory -f make.drivers debug

release:
	@make --no-print-directory -f make.drivers release


clean: clean_debug

clean_debug:
	@make --no-print-directory  -f make.drivers clean_debug

clean_release:
	@make --no-print-directory  -f make.drivers clean_release

clean_all:
	@make --no-print-directory -f make.drivers clean_all

#clean_drivers:
#	@make --no-print-directory -f make.drivers clean_all

