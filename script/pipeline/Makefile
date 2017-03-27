CXXFLAGS =-std=c++0x -ggdb -Wall 
BIN = /rhome/xzhan018/LGT/script/pipeline/bin/
all: ParseHisatSam SingleUnknownAssignCategory
ParseHisatSam: ParseHisatSam.cpp
	g++ $(CXXFLAGS) ParseHisatSam.cpp -o ${BIN}/ParseHisatSam
SingleUnknownAssignCategory: SingleUnknownAssignCategory.cpp
	g++ $(CXXFLAGS) SingleUnknownAssignCategory.cpp -o ${BIN}/SingleUnknownAssignCategory
