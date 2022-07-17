all: Ricker_rb_ratiovar.so Rickerkf.so Rickerkf_fullLL.so

Ricker_rb_ratiovar.so: Ricker_rb_ratiovar.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_rb_ratiovar.cpp')"

Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp','-O1 -g',DLLFLAGS='', framework = 'TMBad')"

Rickerkf_fullLL.so: Rickerkf_fullLL.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf_fullLL.cpp','-O1 -g',DLLFLAGS='', framework = 'TMBad')"

clean:
	rm -rf *o
