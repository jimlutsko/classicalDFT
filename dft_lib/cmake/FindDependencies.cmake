# Find dependencies:

function(find_library_or_set_default env_variable lib_name)

    set (extra_macro_args ${ARGN})

    # Did we get any optional args?
    list(LENGTH extra_macro_args num_extra_args)
    if (${num_extra_args} GREATER 0)
        list(GET extra_macro_args 0 optional_arg)
        set(lib_path ${optional_arg})
        message ("\\--> library searching folder: ${lib_path}")
    endif()

    find_library(${env_variable} NAMES ${lib_name} PATHS ${lib_path} PATH_SUFFIXES "lib")

    if (DEFINED ${env_variable})
        message("\\--> ${lib_name} library found = [" ${${env_variable}} "]")
        set(${env_variable} ${${env_variable}} PARENT_SCOPE)
    else()
        message(WARNING "Lib not found in ${lib_path}: using default ${lib_name}")
        set(${env_variable} ${lib_name} PARENT_SCOPE)
    endif()

    return()
endfunction()
