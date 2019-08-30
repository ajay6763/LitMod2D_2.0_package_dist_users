#----------------------------- tao makefile -----------------------------
#
#Modify options in ./config.mk
#
#Type  'make'  in this directory to compile tao 
#
#tao has been succesfully compiled with this Makefile in: 
# ·iOS 11
# ·IBM AIX Version 3.2 for IBM RISC 6000 workstations.
# ·Linux for pentium processor. 
# ·Hewlett Packard Envizex.
# ·Sun Solaris OS5
#------------------------------------------------------------------------

include config.mk

all:
	(cd src; make)
	@echo; echo; echo Compilation succeeded!
	@(echo "ADD TO YOUR PATH: `pwd`/bin/  AND  `pwd`/script/")
	@(echo "ADD IN .cshrc:    setenv tao_dir `pwd` ")
	@(echo "ADD IN .bashrc:   export tao_dir=`pwd` ")

clean_for_tar:
	(cd src; make clean)
	rm -f src/*.o lib/libreria.o lib/libreria.a 

tao: 
	(cd src; make tao)


vers: 	clean_for_tar
	rm -R -f tao tao_version
	mkdir tao tao/bin
	cp -R -L Makefile config.mk README demo doc include lib script src   tao
	rm -f tao/doc/.first_compilation.txt #tao/lib/sistbanda* version_tmp/lib/surf_proc* version_tmp/lib/thin_sheet*
	tar -chf tao.tar tao
	gzip -f tao.tar
	echo "UPLOADING to github."
	touch tao/bin/touch_something #needed by git add
	mv tao tao_version
	(cd tao_version; git init; git remote add tao https://github.com/danigeos/tao-geo; git add .; git commit -a -mnewVersion; git push -u -f tao master)


upload_version_starting_from_scratch:
	#(git init; git remote add tao https://github.com/danigeos/tao-geo; git add Makefile README bin config.mk demo doc include lib script src; git commit -a -mnewVersion; git push -u -f tao master)


upload:
	#for initialization:  
	#git init; git remote add tao https://github.com/danigeos/tao-geo; git add Makefile README config.mk bin demo doc include lib script src; git rm --cached doc/.first_compilation.txt
	git commit -a 
	git push tao master

