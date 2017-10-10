#output_size
n = 1000

file = open("input.dat", "w", -1)

file.write("{0}\n".format(n))
for i in range(n):
    s = ""
    if (i < n / 5):
        for j in range(n):
            s += "100 "
        s += "\n"
    else:
        for j in range(n):
            if (j > n / 5 * 2 and j < n / 5 * 3):
                s += "100 "
            else:
                s += "0 "

        s += "\n"
    file.write(s)
file.close()