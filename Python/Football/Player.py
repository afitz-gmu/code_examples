#Player.py

class Player (object):
    def __init__ (self, first, last):
        '''the constructor for a player'''
        self.first = first
        self.last = last
        self.rating = 0
        self.info = []


    def update(self, year, team, comp, att, yeards, tds, itcs):
        '''create a list of information for this player for this year and 
append it to the info field. Then call calcrating.'''
        
        self.info.append([year, team, comp, att, yeards, tds, itcs])
        self.calcrating()


    def calcrating(self):
        '''go through all sub-lists in info adding up totals for comps, 
attempts, etc. Then calculate the overall rating for this player. Store it
in the instance variable "rating" '''
        
        comp = 0
        att = 0
        yeards = 0
        tds = 0
        itcs = 0
        
        for line in self.info:
            
            comp += int(line[2])
            att += int(line[3])
            yeards += int(line[4])
            tds += int(line[5])
            itcs += int(line[6])
            
        c = ((((comp / att) * 100) - 30 ) / 20)
        if c > 2.375:
            c = 2.375
                    
        y = ((((yeards / att)) - 3 ) / 4 )
        if y > 2.375:
            y = 2.375
                    
        t = (tds / att) * 20
        if t > 2.375:
            t = 2.375
                
        intcep = (itcs / att) * 25 
        i = 2.375 - intcep  
        
        self.rating = (((c + y + t + i) / 6) * 100)
        
        return self.rating

    def returnName(self):
        '''return the name of the player first last'''
        
        return self.first + " " + self.last


    def returnReverseName(self):
        '''return the name of the player as last, first'''
        
        return self.last + " " + self.first

    def __eq__ (self, other):
        '''determine if this person's name is the same as the other person's
name'''
        
        if self.first == other.first and self.last == other.last:
            return True
        else:
            return False
        
    def __lt__(self,other):
        '''determine if this person's name is less than the other person's
name alphabetically'''
        
        if self.last == other.last:
            if self.first < other.first:
                return True
            else:
                return False
        elif self.last < other.last:
            return True
        else:
            return False

    def __gt__ (self, other):
        '''determine if this person's name is greater than the other person's
name alphabetically'''
        
        if self.last == other.last:
            if self.first < other.first:
                return True
            else:
                return False
        elif self.last > other.last:
            return True
        else:
            return False        
        

    def __str__(self):
        '''return a string of the person's name and their rating in a nice
format'''
        
        playerRating = "{:25}{:>10.2f}".format(self.returnName(), self.rating)
        return playerRating
        
    def calc(self, sublist):
        '''calculate a passer rating for one sub-list year in  the info list'''
        
        comp = int(sublist[2])
        att = int(sublist[3])
        yeards = int(sublist[4])
        tds = int(sublist[5])
        itcs = int(sublist[6])
        
        c = ((((comp / att) * 100) - 30 ) / 20)
        if c > 2.375:
            c = 2.375
                    
        y = ((((yeards / att)) - 3 ) / 4 )
        if y > 2.375:
            y = 2.375
                    
        t = (tds / att) * 20
        if t > 2.375:
            t = 2.375
                
        intcep = (itcs / att) * 25 
        i = 2.375 - intcep
        
        newPassRating = (((c + y + t + i) / 6) * 100)
        return newPassRating
       

    def printInfo(self):
        '''print individual year information about his player including each 
year's passer rating. The list should be in year order. Use calc to assist.'''
        
        self.info.sort()
        print(self.returnName())
        for line in self.info:
            print("{} in {} - {:>6.2f}".format(line[1], line[0], self.calc(line)))
