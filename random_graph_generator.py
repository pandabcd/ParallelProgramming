import random
# To simplify writing in file, we have max_num_edges

s1 = 13
s2 = 12
prob = 0.07
max_num_edges = int(s1*s2*prob)

f = open("random_" + str(s1) + "_" + str(s2) + ".txt", "w")


# f.write(str(s1) + " " + str(s2) + " " + str(max_num_edges) + "\n")

curr_num_edges = 0
for i in range(1,s1+1):
	for j in range(s1+1, s1+s2+1):
		num = random.random()
		if(num<=prob):
			f.write(str(i) + " " + str(j) + "\n")
			curr_num_edges += 1
		if curr_num_edges > max_num_edges:
			exit(0)