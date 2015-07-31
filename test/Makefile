
CPP := /Users/silentvoice/clang+llvm-3.6.2-x86_64-apple-darwin/bin/clang++


ALL : test.cpp
	$(CPP) -g --std=c++14 test.cpp -o test
	./test &> out


.PHONY : clean


clean :
	rm -rf test
	rm -rf test.dSYM
