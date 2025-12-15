all : test_app

test_app : test_app.c
	gcc -o $@ $^ -lm -O3 -g

test : test_app
	./test_app

clean :
	rm -rf *.o test_app

