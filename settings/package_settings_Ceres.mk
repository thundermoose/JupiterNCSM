package_dependencies_Ceres=$(filter-out ./src/Bacchus/programs/%.c .src/Bacchus/settings/%.c,$(shell find .src/Bacchus/ -regex [^\#]*\\.c$)) $(filter-out ./src/Minerva/programs/%.c,$(shell find ./src/Minerva/ -regex [^\#]*\\.c$))

package_compiler_flags_Ceres=-I./src/Bacchus/
