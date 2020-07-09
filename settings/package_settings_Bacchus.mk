test_Bacchus="This is in Bacchus settings"
package_dependencies_Bacchus=$(filter-out ./src/Minerva/programs/%.c,$(shell find ./src/Minerva/ -regex [^\#]*\\.c$))

package_compiler_flags_Bacchus=-I./src/Minerva/
