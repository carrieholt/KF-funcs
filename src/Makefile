all: Ricker_rb_ratiovar.so Rickerkf.so

Ricker_rb_ratiovar.so: Ricker_rb_ratiovar.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_rb_ratiovar.cpp')"

Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp','-O1 -g',DLLFLAGS='', framework = 'TMBad')"

clean:
	rm -rf *o
