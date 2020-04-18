
VERSION = v3

exe_name ?= qfit_$(VERSION)

############################################################################################################
VERSION_DATE  = $(VERSION)_$(shell date "+%d_%m_%Y_%T")
VERSION_DATE  = $(VERSION)_$(shell date "+%d_%m_%Y_%H-%M-%S")
STANDARD = c99
############################################################################################################

GCC ?= g++ 

# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])
# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")
# These flags will override any settings
ifeq ($(i386),1)
	OS_SIZE = 32
	OS_ARCH = i686
endif
ifeq ($(x86_64),1)
	OS_SIZE = 64
	OS_ARCH = x86_64
endif
# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

ifneq ($(DARWIN),) 
      CCFLAGS   := -arch $(OS_ARCH) 
else
  ifeq ($(OS_SIZE),32)
      CCFLAGS   := -O3 -std=c++11
  else
      CCFLAGS   := -O3 -std=c++11
  endif
endif


INC_LIBS =   -lpng -lgsl -lgslcblas `fltk-config --cxxflags --use-forms --use-gl --use-images --ldflags`
#-lGL -lGLU  -lglut -lfltk -lfltk_gl
#INC_LIBS_STATIC = -L/usr/lib/x86_64-linux-gnu/  -static-libgcc  -static-libstdc++ -l:libm.a -l:libpng.a  `fltk-config --cxxflags --ldstaticflags --use-gl --use-images --use-glut`
INC_LIBS_STATIC = -L/usr/lib/x86_64-linux-gnu/ -static-libgcc -static-libstdc++  -l:libgsl.a -l:libgslcblas.a `fltk-config --cxxflags --ldstaticflags --use-images`



#CCFLAGS = $(shell fltk-config --optim)


INC_DIR = $(shell fltk-config  --cxxflags) 

FLTK_LIBS =$(shell fltk-config --use-forms --use-gl --use-images --ldflags)
INC_LIBS = -lpng -lgsl -lgslcblas $(FLTK_LIBS)

FLTK_STATIC_LIBS =$(shell fltk-config --ldstaticflags --use-images)
INC_LIBS_STATIC = -L/usr/lib/x86_64-linux-gnu/ -static-libgcc -static-libstdc++  -l:libgsl.a -l:libgslcblas.a $(FLTK_STATIC_LIBS)





all : $(exe_name)


OBJS := qfit.o ticks.o fit.o plot.o save_image.o jackerr.o
		
MYOBJS  := $(patsubst %.o,%.o,$(OBJS))
deps = $(MYOBJS:.o=.d)

%.o: %.cpp
	@echo "######################### Compiling: "$<" #########################"
	$(VERBOSE)$(GCC) $(CCFLAGS) $(INC_DIR) -MMD -MP   -c $< -o $@ 

$(exe_name): $(OBJS)
	@echo "######################### Building: $(exe_name) #########################"
	$(VERBOSE)$(GCC) $(CCFLAGS)  -o $@ $+ $(INC_LIBS)

static: INC_LIBS = $(INC_LIBS_STATIC)
static: $(exe_name)  

clean:
	rm -f $(OBJS) $(exe_name)

pack: 
	@echo Generating Package qfit_$(VERSION_DATE).tar.gz
	@tar cvfz qfit_$(VERSION_DATE).tar.gz *.cpp *.h  Makefile
	@echo Generated Package qfit_$(VERSION_DATE).tar.gz


-include $(deps)









