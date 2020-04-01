fc = 20
f = open("FC_" + str(fc) + "_" + str(fc) + ".txt", "w")

for i in range(1,fc+1):
	for j in range(1,fc+1):
		f.write(str(i) + " " + str(j+fc) + "\n")