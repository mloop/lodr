arms01a:  arms_main01a.o arms.o
	gcc -o arms01a arms_main01a.o arms.o -lm -lc
arms_main01a.o: arms.h arms_main01a.c
	gcc -c arms_main01a.c -o arms_main01a.o
arms.o: arms.c
	gcc -c arms.c -o arms.o 

