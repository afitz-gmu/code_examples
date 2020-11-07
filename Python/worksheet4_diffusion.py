import matplotlib.pyplot as plt

infile = open('iPad Quarterly Sales 2010 2017.csv', mode = 'r', encoding = None)
header = infile.readline()
data = infile.readlines()
infile.close()

units = []
for line in data:
    line = line.strip()
    date,Number = line.split(',')
    date = str(date)
    Number = float(Number)
    units.append(Number)

lengthTime = len(units)
lengthTime = list(range(0, lengthTime))

plt.plot(lengthTime, units)
plt.xlabel("time")
plt.ylabel("num units (M)")
plt.title("iPad Sales 2010-2017 (M Units)")
plt.show()