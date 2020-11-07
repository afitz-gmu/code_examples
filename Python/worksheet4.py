import matplotlib.pyplot as plt

infile = open('student_data.csv', mode = 'r', encoding = None)
header = infile.readline()
data = infile.readlines()
infile.close()

heightL = []
genderL = []
yearL = []
combinationL = []

for line in data:
    line = line.strip()
    gender,height,year = line.split(',')
    gender = str(gender)
    height = int(height)
    year = str(year)
    heightL.append(height)
    yearL.append(year)
    combination = (gender, height, year)
    combinationL.append(combination)

maleHeight = []
femaleHeight = []
for line in combinationL:
    if line[0] == 'male':
        maleHeight.append(line[1])
    else:
        femaleHeight.append(line[1])

maleSorted = sorted(maleHeight)
lengthM = len(maleSorted)
lengthMale = list(range(0, lengthM))
#plt.plot(lengthMale, maleSorted, label = "male")

femaleSorted = sorted(femaleHeight)
lengthF = len(femaleSorted)
lengthFemale = list(range(0, lengthF))
#plt.plot(lengthFemale, femaleSorted, label = "female")

#plt.xlabel("student number")
#plt.ylabel("height (in)")
#plt.title("Undergraduate Student Heights")

#plt.legend()
#plt.show()

maleFr = []
maleSo = []
maleJr = []
maleSr = []
femaleFr = []
femaleSo = []
femaleJr = []
femaleSr = []
for line in combinationL:
    if line[0] == 'male':
        if line[2] == ' freshman':
            maleFr.append(line[1])
        elif line[2] == ' sophomore':
            maleSo.append(line[1])
        elif line[2] == ' junior':
            maleJr.append(line[1])
        elif line[2] == ' senior':
            maleSr.append(line[1])

    elif line[0] == 'female':
        if line[2] == ' freshman':
            femaleFr.append(line[1])
        elif line[2] == ' sophomore':
            femaleSo.append(line[1])
        elif line[2] == ' junior':
            femaleJr.append(line[1])
        elif line[2] == ' senior':
            femaleSr.append(line[1])

data = [maleFr, maleSo, maleJr, maleSr, femaleFr, femaleSo, femaleJr, femaleSr]
labelL = ['male fr', 'male so', 'male jr', 'male sr', 'female fr', 'female so', 'female jr', 'female sr']
print(data)

plt.boxplot(data, labels = labelL)
plt.ylabel("height (in)")
plt.title("Student Heights by Gender and Class")
plt.show()