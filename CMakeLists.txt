cmake_minimum_required(VERSION 3.16)

project(prophet)

set(CMAKE_CXX_STANDARD 17) # 20 is causing all sorts of warnings
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(/usr/local/include cmdstan/stan/lib/stan_math/lib/tbb_2019_U8/include cmdstan/stan/src cmdstan/lib/rapidjson_1.1.0 cmdstan/lib/CLI11-1.9.1 cmdstan/stan/lib/stan_math cmdstan/stan/lib/stan_math/lib/eigen_3.3.9 cmdstan/stan/lib/stan_math/lib/boost_1.72.0 cmdstan/stan/lib/stan_math/lib/sundials_5.6.1/include)
link_directories(/usr/local/lib)

add_executable(prophet src/prophet.cpp)

target_include_directories(prophet PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(prophet tbb)
target_compile_options(prophet PRIVATE -Wno-unknown-warning-option -Wno-tautological-compare -Wno-sign-compare -D_REENTRANT -Wno-ignored-attributes -O3 -DBOOST_DISABLE_ASSERTS)

find_package(Catch2 REQUIRED)
add_executable(tests tests/src/tests.cpp)
target_include_directories(tests PUBLIC /usr/local/Cellar/catch2/2.13.4/include "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(tests PRIVATE Catch2::Catch2)
