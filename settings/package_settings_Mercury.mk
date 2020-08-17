package_dependencies_Mercury=$(filter-out ./src/Neptune/programs/%.c,$(shell find ./src/Neptune/ -regex [^\#]*\\.c$))
                                                           
package_compiler_flags_Mercury=-I./src/Neptune/
