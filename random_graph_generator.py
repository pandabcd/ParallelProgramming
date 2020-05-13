import random

s1 = 13
s2 = 100
prob = 1

f = open("random_" + str(s1) + "_" + str(s2) + "_" + str(prob) + ".txt", "w")


for i in range(1,s1+1):
	for j in range(s1+1, s1+s2+1):
		num = random.random()
		if(num<=prob):
			f.write(str(i) + " " + str(j) + "\n")