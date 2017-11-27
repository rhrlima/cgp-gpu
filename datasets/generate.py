#Y = X0+X1 + X0*X1 + -X0*(X1*X1) + 0
func = "x*x + x+x + x-x"
#func = "x0+x1 + x0*x1 - x0*x1**2"

in_  = 1
out_ = 1
min_ = 0
max_ = 20

file = open("symbolic2.data", "w")
file.write("{},{},{},\n".format(in_, out_, (max_ - min_)))
for x in range(min_, max_):
	print("{: >5} : {: <10}".format(x, eval(func)))
	file.write("{},{},\n".format(x,  eval(func)))
file.close()
print("Done")