
objects = albite.o anortite.o diopside.o garnet.o Generator.o olivine.o ortenstatitemg90.o spinel.o generator_table_atten_corr.o anelastic_atten_deri_T_P.o
Comp = gfortran
# Flag for compling only. You should not change this
FCLAGS = -c

# Flag for linking. You should not change this
FMLAGS = -o

#### 
Generator_LINUX : $(objects)
	$(Comp) ${FMLAGS} Generator_LINUX $(objects)

%.o : %.for
	$(Comp) ${FCLAGS}  $<
