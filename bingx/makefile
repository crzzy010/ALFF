app : main.o ReadCsvFile.o 
	gcc -o app main.o ReadCsvFile.o

ReadCsvFile.o : readCSVFile.h ReadCsvFile.c 
	gcc -c ReadCsvFile.c 

main.o : main
	gcc -c main.c