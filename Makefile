PROJECT_ROOT_DIR = .
OUTPUT_DIR = bin
OUTPUT_EXE = steroMatching
SOURCE_DIR = $(PROJECT_ROOT_DIR)/src
SOURCE_FILES = $(SOURCE_DIR)/stereo.cpp $(SOURCE_DIR)/FreeImageHelper.cpp $(SOURCE_DIR)/Eigen.h $(SOURCE_DIR)/FreeImageHelper.h
INCLUDE_DIRS = $(PROJECT_ROOT_DIR)/libs/Eigen
LIBRARIES = -lfreeimage
CC_FLAGS = -std=c++14 -w -O3
OUTPUT_PATH = $(OUTPUT_DIR)/$(OUTPUT_EXE)

all:
	mkdir -p $(OUTPUT_DIR)
	g++ $(CC_FLAGS) $(SOURCE_FILES) -o $(OUTPUT_PATH) -I $(INCLUDE_DIRS) $(LIBRARIES) 
	
clean:
	rm -rf $(OUTPUT_DIR)

run:
	./bin/$(OUTPUT_EXE)
