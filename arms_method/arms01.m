arms01:  arms_main01.o arms.o
	gcc -o arms01 arms_main01.o arms.o -lm -lc
arms_main01.o: arms.h arms_main01.c
	gcc -c arms_main01.c -o arms_main01.o
arms.o: arms.c
	gcc -c arms.c -o arms.o 

