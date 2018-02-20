
configure_file(fetch-data.py ${PROJECT_BINARY_DIR}/Data/angiotk_fetch-data.py)

install(PROGRAMS ${PROJECT_BINARY_DIR}/Data/angiotk_fetch-data.py
    DESTINATION bin/)
