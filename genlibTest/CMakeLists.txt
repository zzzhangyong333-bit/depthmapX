set(genlibtest genlibTest)
set(genlibTest_SRCS
    testreadwritehelpers.cpp
    main.cpp
    testsimplematrix.cpp
    testbspnode.cpp
    teststringutils.cpp
    testpafmath.cpp
    testcontainerutils.cpp
    testparallel.cpp)

set(LINK_LIBS
    genlib)

include_directories("../ThirdParty/Catch")

add_executable(${genlibtest} ${genlibTest_SRCS})
find_package(Threads REQUIRED)
target_link_libraries(${genlibtest} ${LINK_LIBS} Threads::Threads)
    
