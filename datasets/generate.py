func = "x*x*x + x*x + x+x"

in_  = 1
out_ = 1
min_ = 0
max_ = 1000000

file_name = "symbolic2_{}.data".format(max_)
file = open(file_name, "w")

file.write("{},{},{},\n".format(in_, out_, (max_ - min_)))

for v in range(min_, max_):
	x = v/max_
	file.write("{},{},\n".format(x/1000,  eval(func)))

file.close()
print("Done")