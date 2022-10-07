include_directories(../)

add_library(gpu_details gpu_details.cpp)
add_library(hybrid hybrid.cpp)
add_library(kernel_manager kernel_manager.cpp)

add_library(accelerator_obj INTERFACE)
target_sources(accelerator_obj INTERFACE $<TARGET_OBJECTS:gpu_details> $<TARGET_OBJECTS:hybrid>
               $<TARGET_OBJECTS:kernel_manager>)

install(TARGETS gpu_details hybrid kernel_manager DESTINATION lib)
install(FILES gpu_details.h hybrid.h hybrid.tpp kernel_manager.h DESTINATION include)
