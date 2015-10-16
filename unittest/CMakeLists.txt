include_directories(gtest/include ../src)
add_executable(unittest scaffoldtest.cpp fibtest.cpp
        nucleotidetest.cpp kmermdtest.cpp kmertest.cpp tstringtest.cpp
        ../src/alignment.cpp ../src/observations.cpp
        ../src/reductions.cpp ../src/coverage.cpp ../src/library.cpp ../src/kmernode.cpp ../src/kmertable.cpp ../src/cliptips.cpp
        ../src/dsnode.cpp ../src/nucleotide.cpp ../src/nodeendstable.cpp
        ../src/settings.cpp ../src/util.cpp ../src/tstring.cpp ../src/kmeroverlap.cpp ../src/graph.cpp
        ../src/suffix_tree.c ../src/solutioncomp.cpp ../src/bubble.cpp ../src/ExpMaxClustering.cpp)

target_link_libraries(unittest readfile fiblib gtest essaMEM
                      gtest_main ${ZLIB_LIBRARIES} ${LAPACK_LIBRARIES}
                      ${GSL_LIBRARIES} pthread)

add_subdirectory(gtest)
