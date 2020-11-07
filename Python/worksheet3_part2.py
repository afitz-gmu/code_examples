import math

class Point:
    """Represents two dimensional point
    Attributes: x, y coordinate"""
    
    def __init__(self, x, y):
        self.x = x 
        self.y = y
        
    def __str__(self):
        print("X coordinate: ", self.x)
        print("Y coordinate: ", self.y)

x = int(input("Please enter an x coordinate: "))
y = int(input("Please enter an y coordinate: "))

point = Point(x, y)

class Circle:
    """creates a circle with a radius and point value"""
    
    def __init__(self, radius):
        self.radius = radius 
        point = Point(x, y)
        
    def print(self):
        print("X coordinate: ", point.x)
        print("Y coordinate: ", point.y)
        print("The radius of the circle is: ", self.radius)
        print("The area of the circle is: %3f" % (self.area))

    def move(self):
        point.x = int(input("Please enter a new x coordinate: "))
        point.y = int(input("Please enter a new y coordinate: "))
        
    def area(self):
        self.area = math.pi * (self.radius * self.radius)
        


radius = float(input("Please enter a radius: "))

circle = Circle(radius)
circle.area()
circle.print()

# calls the move method which will update the x and y values of the circle
circle.move()

#prints out the new circle details
circle.print()
