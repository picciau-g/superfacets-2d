#-------------------------------------------------
# 2014
# Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
# Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
#               Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
#
# Submission to Pacific Graphics 2014
#
#-------------------------------------------------


QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Superfacets
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

QMAKE_MAC_SDK = macosx10.11

macx: {
    QMAKE_MAC_SDK = macosx10.11
}


QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES += main.cpp \
    Vertex3D.cpp \
    Vertex2D.cpp \
    Triangle.cpp \
    Timer.cpp \
    segmenter.cpp \
    Reader.cpp \
    normals.cpp \
    meshvisualizer.cpp \
    Edge.cpp \
    dialogs.cpp \
    colormapper.cpp \
    vertexbasedsegmenter.cpp

HEADERS += \
    Vertex3D.h \
    Vertex2D.h \
    Triangle.h \
    Timer.h \
    Sorting.h \
    segmenter.h \
    Reader.h \
    normals.h \
    meshvisualizer.h \
    Mesh.h \
    Edge.h \
    dialogs.h \
    colormapper.h \
    vertexbasedsegmenter.h \
    common.h


#LIBS += -framework -lGLU
unix:!macx {
LIBS += -lGLU
}
