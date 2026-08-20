if(NOT UNIX)
  return()
endif()

if(
  NOT DEFINED CPACK_TEMPORARY_INSTALL_DIRECTORY
  OR NOT IS_DIRECTORY "${CPACK_TEMPORARY_INSTALL_DIRECTORY}"
)
  message(FATAL_ERROR "CPack did not provide an archive staging directory")
endif()

file(
  GLOB_RECURSE archive_entries
  LIST_DIRECTORIES true
  "${CPACK_TEMPORARY_INSTALL_DIRECTORY}/*"
)
list(PREPEND archive_entries "${CPACK_TEMPORARY_INSTALL_DIRECTORY}")

foreach(archive_entry IN LISTS archive_entries)
  if(IS_SYMLINK "${archive_entry}")
    continue()
  endif()

  if(IS_DIRECTORY "${archive_entry}")
    file(
      CHMOD "${archive_entry}"
      PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
    )
  elseif(archive_entry MATCHES "/bin/AssemblyCpp$")
    file(
      CHMOD "${archive_entry}"
      PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
    )
  else()
    file(
      CHMOD "${archive_entry}"
      PERMISSIONS
        OWNER_READ OWNER_WRITE
        GROUP_READ
        WORLD_READ
    )
  endif()
endforeach()
