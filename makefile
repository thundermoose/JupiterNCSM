.SILENT:
-include makefile.local
ifndef compiler
compiler := gcc -std=gnu99 -Wall -Werror -c
endif

ifndef linker
linker := gcc -std=gnu99 -Wall -Werror
endif

#ifndef thundertester_linker_flags
#thundertester_linker_flags := -lthundertester 
#endif
#
#ifndef thundertester_compiler_flags
#thundertester_compiler_flags := -DTEST
#endif

ifndef test_data_path
test_data_path := "\"./jupiterncsm_test_data/\""
endif

ifndef wigxjpf_path
wigxjpf_path := ../wigxjpf-1.10
endif

ifndef hdf5_comp_flags
hdf5_comp_flags=
endif

ifndef hdf5_link_flags
hdf5_link_flags= -lhdf5
endif

ifndef blas_comp_flags
blas_comp_flags=
endif

ifndef blas_link_flags
blas_link_flags= -lblas -llapack
endif

source_path := src
mode_path := .tmp
object_path := .tmp/objects
dependencies_path := .tmp/dependencies

compiler_flags := -I./$(source_path)/ -I./$(source_path)/Utilities/  -I$(wigxjpf_path)/inc $(hdf5_comp_flags) -fopenmp $(blas_comp_flags)
linker_flags := -lm $(blas_link_flags) -L$(wigxjpf_path)/lib -lwigxjpf $(hdf5_link_flags) -fopenmp -pthread -lconfig 

all_sources := $(shell find ./$(source_path)/ -regex [^\#]*\\.c$)
all_objects := $(all_sources:./$(source_path)/%.c=./$(object_path)/%.o)
dependencies_files := $(all_sources:./$(source_path)/%.c=$(dependencies_path)/%.d)
program_packages_names := $(filter-out Utilities, $(shell ls ./src))
program_package_sources=$(filter ./$(source_path)/$(1)/%.c,$(all_sources))

.SECONDARY: %.x
.PRESIOUS: $(all_objects) $(dependencies_files) %.x

all: $(mode_path)/release.mode $(all_sources)
	make $(program_packages_names:%=release_%)

profiling: $(mode_path)/profiling.mode $(all_sources)
	make $(program_packages_names:%=profiling_%)

debug: $(mode_path)/debug.mode $(all_sources)
	make $(program_packages_names:%=debug_%)

test: $(mode_path)/test.mode $(all_sources)
	make $(program_packages_names:%=test_%)

test_no_logging: $(mode_path)/test_no_logging.mode $(all_sources)
	make $(program_packages_names:%=test_no_logging_%)

clean:
	rm -rf .tmp release debug

$(mode_path)/%.mode:
	rm -rf ./.tmp/
	mkdir -p $(@D)
	echo "Compilation mode $(@F:%.mode=%)" > $@

define link_package
-include settings/package_settings_$1.mk

ifndef package_dependencies_$1
package_dependencies_$1=
endif

ifndef package_compiler_flags_$1
package_compiler_flags_$1=
endif

package_function_sources_$1=$$(filter-out ./$(source_path)/$1/programs/%.c,\
			       $$(filter ./$(source_path)/$1/%.c ./$(source_path)/Utilities/%.c,\
				  $$(all_sources))) $$(package_dependencies_$1)
package_function_objects_$1=$$(package_function_sources_$1:./$(source_path)/%.c=./$(object_path)/%.o) 
package_programs_sources_$1=$$(filter ./$(source_path)/$1/programs/%.c,$$(filter ./$(source_path)/$1/%.c ./$(source_path)/Utilities/%.c,$$(all_sources)))
package_programs_$1=$$(package_programs_sources_$1:./$(source_path)/$1/programs/%.c=$1/%.x)
release_$1: $$(package_programs_$1:%.x=release/%.x)
release/$1/%.x: compiler_flags+=-I./$(source_path)/$1/ $$(package_compiler_flags_$1) -DNDEBUG -DNLOGING -O3 
release/$1/%.x: ./$(object_path)/$1/programs/%.o $$(package_function_objects_$1)
	mkdir -p $$(@D)
	echo "Linking $$@"
	$$(linker) -o $$@ $$^ $$(linker_flags)

profiling_$1: $$(package_programs_$1:%.x=profiling/%.x)
profiling/$1/%.x: compiler_flags+=-I./$(source_path)/$1/ $$(package_compiler_flags_$1) -DNLOGING -pg -ggdb 
profiling/$1/%.x: linker_flags+=-pg
profiling/$1/%.x: ./$(object_path)/$1/programs/%.o $$(package_function_objects_$1)
	mkdir -p $$(@D)
	echo "Linking $$@"
	$$(linker) -o $$@ $$^ $$(linker_flags)

debug_$1: $$(package_programs_$1:%.x=debug/%.x)
debug/$1/%.x: compiler_flags+=-I./$(source_path)/$1/ $$(package_compiler_flags_$1) -ggdb -DDEBUG
debug/$1/%.x: ./$(object_path)/$1/programs/%.o $$(package_function_objects_$1)
	mkdir -p $$(@D)
	echo "Linking $$@ with $$^"
	$$(linker) -o $$@ $$^ $$(linker_flags)

test_$1: $$(package_programs_$1:%.x=test/%.x)
test/$1/%.x: compiler_flags+=-I./$(source_path)/$1/ $$(package_compiler_flags_$1) -ggdb -DDEBUG -DTEST -DTEST_DATA=$(test_data_path)
test/$1/%.x: ./$(object_path)/$1/programs/%.o $$(package_function_objects_$1)
	mkdir -p $$(@D)
	echo "Linking $$@"
	$$(linker) -o $$@ $$^ $$(linker_flags)
#$$@ --run-all-tests

test_no_logging_$1: $$(package_programs_$1:%.x=test_no_logging/%.x)
test_no_logging/$1/%.x: compiler_flags+=-I./$(source_path)/$1/ $$(package_compiler_flags_$1) -ggdb -DDEBUG -DTEST -DTEST_DATA=$(test_data_path) -DNLOGING
test_no_logging/$1/%.x: ./$(object_path)/$1/programs/%.o $$(package_function_objects_$1)
	mkdir -p $$(@D)
	echo "Linking $$@"
	$$(linker) -o $$@ $$^ $$(linker_flags)
endef

$(foreach package, $(program_packages_names), $(eval $(call link_package,$(package))))
# Compilation
-include $(dependencies_files)

$(object_path)/%.o: $(source_path)/%.c $(dependencies_path)/%.d
	echo "Compiling $@ from $<"
	mkdir -p $(@D)
	$(compiler) -o $@ $< $(compiler_flags)

$(dependencies_path)/%.d: $(source_path)/%.c
	mkdir -p $(@D)
	$(compiler) -o /dev/null -MM -MF $@ -MT $(@:$(dependencies_path)/%.d=$(object_path)/%.o) $<
