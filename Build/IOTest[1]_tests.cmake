add_test([=[InTest.BasicTest]=]  [==[C:/Users/andre/Documents/Unief/Master 2/Project/Repo/Build/IOTest.exe]==] [==[--gtest_filter=InTest.BasicTest]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[InTest.BasicTest]=]  PROPERTIES WORKING_DIRECTORY [==[C:/Users/andre/Documents/Unief/Master 2/Project/Repo/Build]==] SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  IOTest_TESTS InTest.BasicTest)
