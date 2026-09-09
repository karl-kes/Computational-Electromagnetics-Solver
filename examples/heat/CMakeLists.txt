set(KESTREL_HEAT_PRECISION "float" CACHE STRING "Heat example precision: float or double")
set_property(CACHE KESTREL_HEAT_PRECISION PROPERTY STRINGS float double)
if(NOT KESTREL_HEAT_PRECISION MATCHES "^(float|double)$")
    message(FATAL_ERROR "KESTREL_HEAT_PRECISION must be float or double")
endif()

set(HEAT_CORE_SOURCES
    src/config/config.cpp
    src/grid/grid.cu
    src/integrator/integrator.cu
    src/simulation/simulation.cu
)
set(HEAT_TEST_SOURCES
    tests/physics.cu
    tests/neumann.cu
    tests/config.cu
    tests/aligned_soa.cu
    tests/initialization.cu
    tests/convergence_order.cu
    tests/convergence.cu
)

if(KESTREL_CUDA)
    # Headers select allocation and inline functions using __CUDACC__.
    # Compile every consumer with the same backend, including .cpp sources.
    set_source_files_properties(${HEAT_CORE_SOURCES} src/main.cpp ${HEAT_TEST_SOURCES}
        PROPERTIES LANGUAGE CUDA
    )
else()
    kestrel_configure_sources(${HEAT_CORE_SOURCES} src/main.cpp ${HEAT_TEST_SOURCES})
endif()

add_library(kestrel-heat-core STATIC ${HEAT_CORE_SOURCES})
target_include_directories(kestrel-heat-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(kestrel-heat-core PUBLIC kestrel::kestrel)
if(KESTREL_HEAT_PRECISION STREQUAL "double")
    target_compile_definitions(kestrel-heat-core PUBLIC HEAT_SOLVER_PRECISION_DOUBLE=1)
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(kestrel-heat-core PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(kestrel-heat-core PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp,-Wall,-Wextra>
)

if(IS_DEBUG)
    target_compile_options(kestrel-heat-core PRIVATE -O0 -g)
    if(HAS_SANITIZERS)
        target_compile_options(kestrel-heat-core PUBLIC
            -fsanitize=address,undefined -fno-omit-frame-pointer
        )
        target_link_options(kestrel-heat-core PUBLIC -fsanitize=address,undefined)
    endif()
else()
    target_compile_options(kestrel-heat-core PRIVATE -O3)
endif()

if(KESTREL_CUDA)
    set_target_properties(kestrel-heat-core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
elseif(ipo_supported AND NOT IS_DEBUG)
    set_target_properties(kestrel-heat-core PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_executable(kestrel-heat src/main.cpp)
target_link_libraries(kestrel-heat PRIVATE kestrel-heat-core)

if(BUILD_TESTS)
    foreach(test_source IN LISTS HEAT_TEST_SOURCES)
        get_filename_component(test_name "${test_source}" NAME_WE)
        add_executable(kestrel-heat-test-${test_name} "${test_source}")
        target_link_libraries(kestrel-heat-test-${test_name} PRIVATE kestrel-heat-core)
        set(test_args)
        if(test_name STREQUAL "convergence")
            set(test_args --nx 16 --ny 16 --nz 16 --steps 2 --output-interval 0)
        endif()
        add_test(NAME heat-${test_name} COMMAND kestrel-heat-test-${test_name} ${test_args})
        set_tests_properties(heat-${test_name} PROPERTIES
            TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=4" LABELS heat
        )
        if(KESTREL_CUDA)
            set_tests_properties(heat-${test_name} PROPERTIES SKIP_RETURN_CODE 77)
        endif()
    endforeach()

    add_test(NAME heat-invalid-nx COMMAND kestrel-heat --nx 2 --steps 0)
    add_test(NAME heat-invalid-alpha COMMAND kestrel-heat --alpha 0 --steps 0)
    add_test(NAME heat-invalid-dx COMMAND kestrel-heat --dx 0 --steps 0)
    set_tests_properties(heat-invalid-nx heat-invalid-alpha heat-invalid-dx
        PROPERTIES WILL_FAIL TRUE TIMEOUT 10 LABELS heat
    )
endif()
