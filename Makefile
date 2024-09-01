all:
	./genRandomChain
	gcc -o genRandomChain genRandomChain.c -lm
	./genRandomChain 50 1 1.53 0.25 50 950 50 650 370 1600 110 10
	gcc -o genCrystalAndChains genCrystalAndChains.c -lm
	./genCrystalAndChains -12 12 -12 12 -6 6 -100 100 -100 100 -100 100 10 10 10 10 1.53 100 1
