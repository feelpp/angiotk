
configure_file(fetch-data.py ${PROJECT_BINARY_DIR}/Data/fetch-data.py)

install(PROGRAMS ${PROJECT_BINARY_DIR}/Data/fetch-data.py
    DESTINATION bin/AngioTK/Data/)
