# fc = full connectivity
# Format - ith index of First row tells about the index of first edge i.
fc = 2
f = open("FC_" + str(fc) + "_" + str(fc) + ".txt", "w")

for i in range(0, fc*fc, fc):
	f.write(str(i) + " ")

f.write("\n")

for i in range(1,fc+1):
	for j in range(1,fc+1):
		f.write(str(i) + " " + str(j+fc) + "\n")