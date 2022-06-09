all: Ricker_rb_ratiovar.so Rickerkf.dll

Ricker_rb_ratiovar.so: Ricker_rb_ratiovar.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_rb_ratiovar.cpp','-O1 -g', DLLFLAGS='', framework = 'TMBad')"

Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp', '-O1 -g',DLLFLAGS='', framework = 'TMBad')"

clean:
	rm -rf *.so
