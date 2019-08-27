arms02:  arms_main02.o arms.o
	gcc -o arms02 arms_main02.o arms.o -lm -lc
arms_main02.o: arms.h arms_main02.c
	gcc -c arms_main02.c -o arms_main02.o
arms.o: arms.c
	gcc -c arms.c -o arms.o 

