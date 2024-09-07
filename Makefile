all:
	./genRandomChain
	gcc -o genRandomChain genRandomChain.c -lm
	./genRandomChain 50 1 1.53 0.25 50 950 50 650 370 1600 110 10
	gcc -o genCrystalAndChains genCrystalAndChains.c -lm
	./genCrystalAndChains -10 10 -30 30 -30 30 -100 100 -100 100 -100 100 10 10 10 30 1.53 100 1
	./genCrystalAndChains -10 10 -15 15 -15 15 -100 100 -100 100 -100 100 10 10 10 10 1.53 100 1
