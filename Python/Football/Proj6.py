#Proj6.py
#quarterback rating using classes
#
#
#This program take a csv file of quarterbacks and creates a player dictionary.
#It then uses the Player class to calculate the overall passer rating which is
#updated in the player dictionary.  The input allows for the user to enter
#a specific quarterback's name.  Then, if he exists it will ask whether or not
#the user wants the overall career passer rating or the passer rating by each 
#season chronologically which is returned by the player class.

from Player import Player

infoList = []
playerDict = {}

        
#-------------------------------------------------------------------------------

def createDict():
    '''Runs through the csv file, strips away all punctuation and appends
    the lines of the csv file to a dictionary that is created using the player 
    class.  This dictionary contains the players names and their overall career
    rating'''    
    
    file = open("passers.csv", "r")
    file.readline()    

    for line in file:
        
        if line != "":
            
            line = line.strip().split(",")
            playerName = line[0] + " " + line[1]
            
            if playerName not in playerDict.keys():
                playerDict[playerName] = Player(line[0], line[1])
            playerDict[playerName].update(int(line[4]), line[3], float(line[6]),\
                float(line[7]), float(line[8]), float(line[9]), float(line[12]))   
            
    file.close()    
        
    return playerDict


#-------------------------------------------------------------------------------


playerDict = createDict()
                                  
#Creates a list of all of the players in the playerDict and their overall career
#passer rating
for line in playerDict.values():
    infoList.append(line)

#Sorts the list by lastname and prints the list
for line in sorted(infoList):
    print(line)
    
print()
ask = input("Do you want to information about a player? ")
print(ask)
ask = ask.lower()
print()

while ask == "y":
    
    playerName = input("Enter a player's name: ")
    print(playerName)
    print()
    
    #If the QB entered is in the dictionary then it asks what information the 
    #user wants.
    if playerName in playerDict.keys():
        
        print("Pick one")
        print("a) Overall passer ratings")
        print("b) Individual years and ratings")
        print()
        
        userChoice = input("Enter choice: ")
        print(userChoice)
        userChoice = userChoice.lower()
        
        if userChoice == "a":
            
            print()
            print(playerDict[playerName])
        
        elif userChoice == "b":
            
            print()
            playerDict[playerName].printInfo()
        
        else:
            
            print()
            print("You've entered an illegal choice")
        
    else:
        
        print("This player is not in the system.")
    
    #Continues to prompt until the user no longer wants to search for another QB
    print()
    ask = input("Are you interested in another player? ")
    print(ask)
    print()
    ask = ask.lower()
    

