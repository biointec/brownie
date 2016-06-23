add_executable(brownie kmercounttable.cpp  kmeroverlaptable.cpp readcorrection.cpp alignment.cpp coverage.cpp library.cpp kmernpp.cpp kmertable.cpp correctgraph.cpp dsnode.cpp nucleotide.cpp nodeendstable.cpp settings.cpp util.cpp tstring.cpp kmeroverlap.cpp graph.cpp brownie.cpp refcomp.cpp suffix_tree.c)

target_link_libraries(brownie readfile essaMEM pthread)

if (ZLIB_FOUND)
   target_link_libraries(brownie ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS brownie RUNTIME DESTINATION bin)
add_subdirectory(readfile)
add_subdirectory(essaMEM-master)
