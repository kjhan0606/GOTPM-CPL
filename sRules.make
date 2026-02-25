#CAMBTYPE = CAMB-1.5.3
CAMBTYPE = CAMB

ifeq ($(CAMBTYPE), CAMB-1.5.3)
	CAMBDIR = ./CAMB-1.5.3/fortran/Release
	CAMBDIR2 = ./CAMB-1.5.3/forutils/Release
	CAMBDIRLIB = ./CAMB-1.5.3/fortran/Release
	CAMBDIRLIB2 = ./CAMB-1.5.3/forutils/Release
	SUBDIRS =
	MYSUB = mysub.1.5.3
	VIEWPOWER = viewpower.1.5.3
	MKPOWER = mkpower.1.5.3
	# This is VAR_DE
	CAMBDE = 1
	CAMBL2 = -lforutils
	# This is for the input CDM P(k)
#	CAMBOPTION = -D_CDM_PK_
 	CAMBOPTION = -D_TOTM_PK_
else
	CAMBDIR = ./CAMB
	CAMBDIR2 = ./CAMB
	CAMBDIRLIB = ./CAMB
	CAMBDIRLIB2 = ./CAMB
	SUBDIRS = $(CAMBDIR)
	MYSUB = mysub
	MKPOWER = mkpower
	VIEWPOWER = viewpower
	# This is LCDM
	CAMBDE = 0
	CAMBL2 =
	# This is for the input total Matter P(k)
	CAMBOPTION = -D_TOTM_PK_
endif

CAMBOPTION += -DDEMODEL=$(CAMBDE) 
CAMBL = -lcamb  $(CAMBL2)
