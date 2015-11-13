include_directories(gtest/include ../src)
add_executable(unittest alignmenttest.cpp scaffoldtest.cpp fibtest.cpp readfiletest.cpp
        nucleotidetest.cpp kmermdtest.cpp kmertest.cpp tstringtest.cpp
        ../src/tstring.cpp ../src/nucleotide.cpp ../src/kmeroverlap.cpp ../src/alignment.cpp)

target_link_libraries(unittest readfile fiblib gtest essaMEM
                      gtest_main ${ZLIB_LIBRARIES} ${LAPACK_LIBRARIES}
                      ${GSL_LIBRARIES} pthread)

add_subdirectory(gtest)
