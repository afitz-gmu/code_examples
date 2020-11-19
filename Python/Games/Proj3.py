#
#This program contains five games. The first game finds all words at a 
#particular length, with only one vowel, and excluding a certain letter. The 
#second game finds all words in the dictionary that contain all the letters of
#a word inputted by the user. The third finds all words at a particular length
#that consist of only consecutive state abbreviations. The fourth finds all 
#words that contain a string inputted by the user. The fifth finds all 
#palindromes at a particular length.

dictionary = open("dictionary.txt", "r")


def exclude_char(word):
    """Return words that are a certain length, exclude a certain character 
    and have one vowel."""
    
    length = input("Please enter the word length you are looking for: ")
    print(length)
    wordLength = int(length)
    exclude = input("Please enter the letter you'd like to exclude: ")
    print(exclude) 
    print()
    word_count = 0
    
    for word in dictionary:
        vowel = "aeiou"
        count = 0
        count_vowel = 0
        word = word.strip()
        
        for char in word:
            
            if char in vowel:
                count_vowel = count_vowel + 1
            
            if char in exclude:
                    count = count + 1
        
        if len(word) == wordLength and count_vowel <= 1 and count == 0:
            print(word)
            word_count  += 1
    
    if word_count == 0:
        print("There are no words that fit this criteria.")

#-------------------------------------------------------------------------------
    
def same_word(word):
    """Returns all words that are under a maximum length and contain all of
    the characters of a certain word that was entered""" 
    
    word_selected = input("Enter word: ")
    print(word_selected)
    length = input("What is the maximum length of the words you want: ")
    print(length)
    print()
    wordLength = int(length)
    word_count = 0 
    
    for word in dictionary:
        temp = word_selected
        word = word.strip()
        
        #Sees if the length of the word is less than or equal to the max length
        if len(word) <= wordLength:
            
            for char in word:
                
                #Checks to see if that character is also in temp
                if char in temp:
                    temp = temp.replace(char,"",1)
        
        if temp == "":
            print(word)
            word_count += 1
    
    if word_count == 0:
        print("There are no words that fit this criteria.")    

#-------------------------------------------------------------------------------       

def state_abreviations(word):
    """Returns all words that are made up of overlapping state abreviations
    and are a certain length"""    
    states = ['AL','AK','AZ','AR','CA','CO','CT','DE','DC','FL','GA','HI',\
'ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT',\
'NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD',\
'TN','TX','UT','VT','VA','WA','WV','WI','WY']
    length = input("Enter the length of the words you are looking for: ")
    print(length)
    print()
    word_length = int(length)
    word_count = 0
    
    for word in dictionary:
        count_one = 0 
        count_two = 2
        word = word.strip().upper()
       
        if len(word) == word_length:
            abreviation = word[count_one:count_two] 
            count_one += 1
            count_two += 1            
           
            #Checks the pairs of characters to make sure they are state abrev.
            while abreviation in states:
                abreviation = word[count_one:count_two] 
                count_one += 1
                count_two += 1                    
                
                #If the count one gets to the end the whole word is state abrev.
                if count_one == word_length:
                    print(word)
                    word_count += 1
    
    if word_count == 0:
        print("There are no words that fit this criteria.")     

#-------------------------------------------------------------------------------

def same_string(word):
    """Returns all words that contain the string we are looking for"""    
    phrase = input("Please enter the string you are looking for: ")
    print(phrase)
    print()
    word_count = 0
    
    for word in dictionary:
        word = word.strip()        
            
        #Checks if the phrase entered is in any words in the dictionary
        if phrase in word:
            print(word)
            word_count += 1
                
    
    if word_count == 0:
        print("There are no words that fit this criteria.")     
        
#-------------------------------------------------------------------------------

def palindrome(word):
    """Returns all palindromes that are a certain length"""
    size = input("Enter the length of the palindromes you desire: ")
    print(size)
    print()
    palindrome_length = int(size)
    word_count = 0    
    
    for word in dictionary:
        word = word.strip()
        
        #Checks the word in the dictionary against itself reversed
        if len(word) == palindrome_length and word == word[::-1]:
            print(word)
            word_count += 1
    
    if word_count == 0:
        print("There are no words that fit this criteria.")     

#-------------------------------------------------------------------------------

#Main Program
choice = ""
print()
print("Lets play a game")
while choice != "q":
    dictionary.seek(0)
    print()
    print("Choose which game you want to play")
    print("a) Find words with only one vowel and excluding a specific letter.")
    print("b) Find words containing all the letters of another word with a maximum length")
    print("c) Find words containing overlapping state abbreviations")
    print("d) Find words containing an exact string of characters")
    print("e) Find all the palindromes of a particular length")
    print("q) quit")
    print()
    
    choice = input("Enter a choice: ")
    print(choice)
    print()
    choice = choice.lower()
    
    if choice == 'a':
        exclude_char(dictionary)
        
    elif choice == 'b':
        same_word(dictionary)
        
    elif choice == 'c':
        state_abreviations(dictionary)
        
    elif choice == 'd':
        same_string(dictionary)
        
    elif choice == 'e':
        palindrome(dictionary)
    
    elif choice == 'q':
        print("Thanks for playing")  
        dictionary.close()
        
    elif choice != 'a' or 'b' or 'c' or 'd' or 'e' or 'q':
        print("You've entered an incorrect choice. Try again")
        
    
