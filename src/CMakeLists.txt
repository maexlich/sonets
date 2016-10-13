include_directories(${GSL_INCLUDE_DIR})
SET(GSL_LIBS ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

file(GLOB LIB_SRC
        "lib/*.hpp"
        "lib/*.cpp"
        )
ADD_LIBRARY(LIB OBJECT ${LIB_SRC})

add_executable(run_secorder run_secorder.cpp $<TARGET_OBJECTS:LIB>)
target_link_libraries(run_secorder ${GSL_LIBS})

add_executable(run_secorder_2p run_secorder_2p.cpp $<TARGET_OBJECTS:LIB>)
target_link_libraries(run_secorder_2p ${GSL_LIBS})