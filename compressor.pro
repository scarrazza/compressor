TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += src/main.cc \
    src/utils.cc \
    src/Minimizer.cc \
    src/Estimators.cc \
    src/Grid.cc

HEADERS += inc/utils.hh \
    inc/Minimizer.hh \
    inc/Estimators.hh \
    inc/Grid.hh \

INCLUDEPATH = inc/ $$system(lhapdf-config --incdir) $$system(root-config --incdir)
LIBS += $$system(lhapdf-config --libs) $$system(gsl-config --libs) $$system(root-config --glibs)

OBJECTS_DIR = obj
