import math
import matplotlib.pyplot as plt


class ZeroVectorError(Exception):
    "Tried to perform operation on zero vector"
    pass

class Vector: # 2-vector over R

    def __init__(self,xcomp,ycomp): # constructor based on x and y components
        self._vector = [None,None] # Initialising attributes

        self.setXComp(xcomp)
        self.setYComp(ycomp)

    def __str__(self):
        return "X component {} : Y component {}".format(self.getXComp(), self.getYComp())

    def setXComp(self,xcomp):
        if not isinstance(xcomp, (float,int)):
            raise TypeError
        else:
            self._vector[0] = float(xcomp)

    def setYComp(self,ycomp):
        if not isinstance(ycomp, (float,int)):
            raise TypeError
        else:
            self._vector[1] = float(ycomp)

    def getXComp(self) -> float:
        return self._vector[0]

    def getYComp(self) -> float:
        return self._vector[1]

    def getTuple(self) -> tuple: # returns x and y components as tuple
        return (self.getXComp(),self.getYComp())

    def getMagnitude(self) -> float:
        x = self.getXComp()
        y = self.getYComp()
        if x == 0 and y == 0:
            return 0
        else:
            return math.sqrt(x ** 2 + y ** 2)

    def getDirection(self) -> float: # returns angle between positive x axis and direction, in Radians. Based upon arguments of complex numbers
        x = self.getXComp()
        y = self.getYComp()
        if x != 0:
            fraction = abs(y/x)

        if x > 0:
            if y >= 0: #first quadrant
                direction = math.atan(fraction)
            else: # fourth quadrant
                direction = -1 * math.atan(fraction)

        elif x < 0:
            if y>= 0: #second quadrant
                direction = math.pi - math.atan(fraction)
            else: # third quadrant
                direction = -1 * (math.pi - math.atan(fraction))

        else: # x == 0
            if y > 0:
                direction = math.pi/2 # straight up
            elif y < 0:
                direction = -1 * math.pi/2 # straight down
            else:
                raise ZeroVectorError # Zero Vectors have no direction

        return direction

"""def vectorTest(x,y): # Test subroutine
v = Vector(x,y)
print(v.getXComp())
print(v.getYComp())
print(v.getMagnitude())
print(v.getDirection())"""


class Projectile:

    def __init__(self, v_x, v_y, mass, area, dragCoefficient):
    #Initialising attributes
        self._velocity = Vector(v_x,v_y) #m/s
        self._mass = None #kg
        self._area = None #m^2
        self._Cd = None #dimensionless
        self._displacement = Vector(0,0) # Projectile starts from origin, m

        self.setMass(mass)
        self.setArea(area)
        self.setCd(dragCoefficient)

    def __str__(self):
        return "Velocity {} : Mass {} : Area {} : Drag Coefficient {} : Displacement {}".format(self.getVelocity(), self.getMass(), self.getArea(), self.getCd(), self.getDisplacement())

    def setVelocity(self, velocity):
        if not isinstance(velocity, Vector):
            raise TypeError
        else:
            self._velocity = velocity

    def setMass(self, mass):
        if not isinstance(mass, (int,float)):
            raise TypeError
        elif mass == 0:
            raise ValueError
        else:
            self._mass = float(mass)

    def setArea(self, area):
        if not isinstance(area, (int,float)):
            raise TypeError
        elif area == 0:
            raise ValueError
        else:
            self._area = float(area)

    def setCd(self, dragCoefficient):
        if not isinstance(dragCoefficient, (int,float)):
            raise TypeError
        else:
            self._Cd = float(dragCoefficient)

    def setDisplacement(self, displacement):
        if not isinstance(displacement, Vector):
            raise TypeError
        else:
            self._displacement = displacement

    def getVelocity(self) -> Vector:
        return self._velocity

    def getMass(self) -> float:
        return self._mass

    def getArea(self) -> float:
        return self._area

    def getCd(self) -> float:
        return self._Cd

    def getDisplacement(self) -> Vector:
        return self._displacement

"""def projectileTest(): # Test for projectile
velocity = Vector(0,-1)
displacement = Vector(4,3)
mass = 5
area = 3
dragCoefficient = 6
projectile = Projectile(velocity, mass, area, dragCoefficient)
print(projectile.getVelocity())
print(projectile.getMass())
print(projectile.getArea())
print(projectile.getCd())
projectile.setDisplacement(displacement)
print(projectile)

projectileTest()"""

class Environment:

    def __init__(self,mass,airDensity,radius):
    # Initialising values
        self._gravConst = 6.67430 * 10**-11 #m^3 kg^-1 s^-2 (National Institute of Standards and Technology, 2019)
        self._mass = None # kg
        self._airDensity = None # kg/m^3
        self._radius = None # m

        self.setMass(mass)
        self.setAirDensity(airDensity)
        self.setRadius(radius)

    def __str__(self):
        return "Gravitational Constant {} : Mass {} : Air Density {} : Radius {}".format(self.getGravConst(), self.getMass(), self.getAirDensity(), self.getRadius())

    def setMass(self, mass):
        if not isinstance(mass, (int,float)):
            raise TypeError
        elif mass == 0:
            raise ValueError
        else:
            self._mass = float(mass)

    def setAirDensity(self, airDensity):
        if not isinstance(airDensity, (int,float)):
            raise TypeError
        
        else:
            self._airDensity = float(airDensity)

    def setRadius(self, radius):
        if not isinstance(radius, (int,float)):
            raise TypeError
        elif radius == 0:
            raise ValueError
        else:
            self._radius = float(radius)

    def getGravConst(self) -> float:
        return self._gravConst

    def getMass(self) -> float:
        return self._mass

    def getAirDensity(self) -> float:
        return self._airDensity

    def getRadius(self) -> float:
        return self._radius

class Simulator:

    def __init__(self, timeStep, speed, projectileMass, area, dragCoefficient, environmentMass, airDensity, radius):
        self._timeStep = None # s quantised time unit to approximate continuous time
        self._speed = None # m/s
        self._mass = None #kg
        self._area = None #m^2
        self._Cd = None #dimensionless

        self.setTimeStep(timeStep)
        self.setSpeed(speed)
        self.setMass(projectileMass)
        self.setArea(area)
        self.setCd(dragCoefficient)
        self._environment = Environment(environmentMass,airDensity,radius)

    def __str__(self):
        return "Time Step {} : Speed {} : Environment {}".format(self.getTimeStep(), self.getSpeed(), self.getEnvironment())

    def setAngle(self, angle): # for trial
        if not isinstance(angle, (float, int)):
            raise TypeError
        elif not(0 <= angle <= (math.pi)/2):
            raise ValueError
        else:
            return float(angle)

    def setTimeStep(self, timeStep):
        if not isinstance(timeStep, (int,float)):
            raise TypeError
        elif timeStep == 0:
            raise ValueError
        else:
            self._timeStep = float(timeStep)

    def setSpeed(self, speed):
        if not isinstance(speed, (float,int)):
            raise TypeError
        elif speed < 0:
            raise ValueError
        else:
            self._speed = float(speed)

    def setEnvironment(self, environment):
        if not isinstance(environment, Environment):
            raise TypeError
        else:
            self._environment = environment

    def getTimeStep(self) -> float:
        return self._timeStep

    def getSpeed(self) -> float:
        return self._speed

    def getEnvironment(self) -> Environment:
        return self._environment

    #From Projectile Class, to prevent repeatedly passing arguments into launch

    def setMass(self, mass):
        if not isinstance(mass, (int,float)):
            raise TypeError
        elif mass == 0:
            raise ValueError
        else:
            self._mass = float(mass)

    def setArea(self, area):
        if not isinstance(area, (int,float)):
            raise TypeError
        elif area == 0:
            raise ValueError
        else:
            self._area = float(area)

    def setCd(self, dragCoefficient):
        if not isinstance(dragCoefficient, (int,float)):
            raise TypeError
        else:
            self._Cd = float(dragCoefficient)

    def getMass(self) -> float:
        return self._mass

    def getArea(self) -> float:
        return self._area

    def getCd(self) -> float:
        return self._Cd

    def launch(self, angle, fullListReturn = True) -> list: # produces list of coordinates during flight, angle in Radians. For options 2 and 3, only the last data value is required, so only the last value is appended.
        time = self.getTimeStep()
        displacements = [(0,0)]
        previousDisplacement = (0,0)
        # instantiate and get references to objects
        projectile = Projectile(self.getSpeed() * math.cos(angle), self.getSpeed() * math.sin(angle), self.getMass(), self.getArea(), self.getCd())
        environment = self.getEnvironment()
        displacement = projectile.getDisplacement()
        velocity = projectile.getVelocity()
        #Initial movement before forces exerted, ensuring that projectiles fired at 0 degrees have 0 displacement
        if angle != 0:
            displacement.setXComp(displacement.getXComp() + velocity.getXComp()*time)
            displacement.setYComp(displacement.getYComp() + velocity.getYComp()*time)
        while displacement.getYComp() > 0: #travels at velocity for the time in timeStep, then recalculates forces and velocity. Repeats until is at or below y=0
            if fullListReturn:
                displacements.append(displacement.getTuple())
            
            dragMagnitude = 0.5 * environment.getAirDensity() * projectile.getArea() * projectile.getCd() * (velocity.getMagnitude() ** 2) # D = .5 * rho * A * Cd * v^2

            weight = (environment.getGravConst() * environment.getMass() * projectile.getMass())/((environment.getRadius() + displacement.getYComp())**2) # F = GMm/r^2

            dragAngle = velocity.getDirection() + math.pi # opposes motion, opposite direction therefore pi radians to projectile direction
            dragXComp = dragMagnitude * math.cos(dragAngle)
            dragYComp = dragMagnitude * math.sin(dragAngle)

            v_x = velocity.getXComp() + (dragXComp * time)/projectile.getMass() # v = u + (Ft)/m
            v_y = velocity.getYComp() + (dragYComp - weight)*time/projectile.getMass()
            velocity.setXComp(v_x)
            velocity.setYComp(v_y)

            previousDisplacement = displacement.getTuple()
            displacement.setXComp(displacement.getXComp() + velocity.getXComp()*time) # s = vt
            displacement.setYComp(displacement.getYComp() + velocity.getYComp()*time)

        if displacement.getYComp() == 0:
            displacements.append(displacement.getTuple())
        else:
            proportion = displacement.getYComp() / (velocity.getYComp() * time) # Proportion of final movement above ground
            s_x = previousDisplacement[0] + proportion * velocity.getXComp() * time # Intersection with ground
            displacements.append([s_x,0])
        print("Launch Complete")
        return displacements

    def FindOptimalAngle(self, fullListReturn = True) -> tuple: # returns list of angles and corresponding displacements, the greatest displacement and the greatest angle. For option 3, only the greatest angle is required, so the list is left empty.
        lowerBound = 0
        upperBound = math.pi / 2
        returnList = []
        greatestDisplacement, greatestAngle = 0,0
        while upperBound-lowerBound >= 0.00001: # absolute error, to prevent long run times
            increment = (upperBound - lowerBound)/10
            displacementsList = [self.launch(angle, False) for angle in [lowerBound + increment * i for i in range(11)]]
            for i in range(0,11):
                displacement,angle = displacementsList[i][-1][0],lowerBound + increment * i
                if fullListReturn:
                    returnList.append([angle,displacement]) # in x-y format
                if displacement > greatestDisplacement:
                    greatestDisplacement, greatestAngle = displacement, angle
            lowerBound, upperBound = greatestAngle - increment, greatestAngle + increment
        if fullListReturn:
            returnList.sort() # sorts list into ascending order by angle
        return (returnList, greatestDisplacement, greatestAngle)


def checkFloat(trial) -> float: # ensures entered value is float, raises error immediately if incorrect value entered
    if not isinstance(trial, (int,float,str)):
        raise TypeError
    else:
        try:
            trial = float(trial)
        except:
            raise TypeError
        else:
            return float(trial)

def inputSciNotation(mantissa,power) -> float: # allows user to enter scientific notation
    return float(mantissa * 10 ** power)

def main():
    print("Welcome to the projectile simulator")
    option = 0
    while option != 4:
        print("What would you like to do? Please enter the corresponding number")
        print("1: Simulate single launch for given conditions")
        print("2: Find optimal angle for given conditions")
        print("3: Investigate change in optimal angle for varying conditions")
        print("4: Exit")
        option = int(input())
        if option == 1:
            option1()
        elif option == 2:
            option2()
        elif option == 3:
            option3()
        elif not (isinstance(option, int) and 1 <= option <= 4):
            print("That is not a valid input")
        else:
            print("End of program")

def inputValues(variable = None, value = None) -> tuple: # will use variable and value from option 3 if given, returns list of parameters for simulator object instantiation
    print("Enter the mantissa as a number, using SI units. Then enter the power of 10.")
    if variable != 1:
        timeStep = inputSciNotation(checkFloat(input("Enter time increment the simulator will use, in s \n")),checkFloat(input()))
    else:
        timeStep = value
    if variable != 2:
        speed = inputSciNotation(checkFloat(input("Enter initial projectile speed, in m/s \n")),checkFloat(input()))
    else:
        speed = value
    if variable != 3:
        projectileMass = inputSciNotation(checkFloat(input("Enter projectile mass, in kg \n")),checkFloat(input()))
    else:
        projectileMass = value
    if variable != 4:
        area = inputSciNotation(checkFloat(input("Enter projectile cross-sectional area, in m^2 \n")),checkFloat(input()))
    else:
        area = value
    if variable != 5:
        Cd = inputSciNotation(checkFloat(input("Enter projectile drag coefficient \n")),checkFloat(input()))
    else:
        Cd = value
    if variable != 6:
        environmentMass = inputSciNotation(checkFloat(input("Enter planet mass, in kg \n")),checkFloat(input()))
    else:
        environmentMass = value
    if variable != 7:
        airDensity = inputSciNotation(checkFloat(input("Enter fluid density, in kg/m^3 \n")),checkFloat(input()))
    else:
        airDensity = value
    if variable != 8:
        radius = inputSciNotation(checkFloat(input("Enter planet radius, in m \n")),checkFloat(input()))
    else:
        radius = value
    return timeStep, speed, projectileMass, area, Cd, environmentMass, airDensity, radius

def option1():
    simulator = Simulator(*inputValues()) # *-operator unpacks lists and tuples
    angle = checkFloat(input("Enter the angle of launch in degrees \n"))
    displacements = simulator.launch(angle * math.pi / 180) # converts into radians
    plt.xlabel("Horizontal Displacement / m")
    plt.ylabel("Vertical Displacement / m")
    plt.title("Projectile Trajectory Under Given Conditions")
    plt.grid(True)
    xValues = [pair[0] for pair in displacements]
    yValues = [pair[1] for pair in displacements]
    plt.plot(xValues,yValues)
    plt.show()
    print("The range was {} metres".format(displacements[-1][0]))

def option2():
    simulator = Simulator(*inputValues())
    values, greatestDisplacement, greatestAngle = simulator.FindOptimalAngle()
    plt.xlabel("Angle / Degrees")
    plt.ylabel("Horizontal Displacement / m")
    plt.title("Relationship Between Angle of Launch and Horizontal Displacement Under Given Conditions")
    plt.grid(True)
    xValues = [pair[0] * (180/math.pi) for pair in values]
    yValues = [pair[1] for pair in values]
    plt.plot(xValues,yValues)
    plt.show()
    print("The greatest displacement was {} metres, when projected at an angle of {} degrees".format(greatestDisplacement, greatestAngle * (180/math.pi)))

def option3():
    print("Enter the number corresponding to the value you would like to vary:")
    print("1: Simulator time interval")
    print("2: Initial projectile speed")
    print("3: Projectile mass")
    print("4: Projectile cross-sectional area")
    print("5: Projectile drag coefficient")
    print("6: Planet mass")
    print("7: Atmospheric fluid density")
    print("8: Planet radius")
    variable = int(input())
    print("Enter the lower bound's mantissa, in SI units. Then enter the power of 10. Repeat for the upper bound")
    lowerBound = inputSciNotation(checkFloat(input("Lower Bound: ")),checkFloat(input()))
    upperBound = inputSciNotation(checkFloat(input("Upper Bound: ")),checkFloat(input()))
    simulator = Simulator(*inputValues(variable, lowerBound))
    increment = (upperBound - lowerBound)/20 # Increment for range
    xValues = [lowerBound + i * increment for i in range(21)] # includes lower and upper bound
    yValues = []
    functionMap = {1: simulator.setTimeStep,
    2: simulator.setSpeed,
    3: simulator.setMass,
    4: simulator.setArea,
    5: simulator.setCd,
    6: simulator.getEnvironment().setMass,
    7: simulator.getEnvironment().setAirDensity,
    8: simulator.getEnvironment().setRadius} # dictionary mapping input to function call

    i = 1
    for value in xValues:
        functionMap.get(variable)(value) # calls function corresponding to variable passing the value as an argument
        yValues.append(simulator.FindOptimalAngle(False)[2] * 180/math.pi) # only appends optimal angle, not other returned values
        print("Iteration {} complete".format(i))
        i+=1


    stringMap = {1: "Simulator time interval",
    2: "Initial projectile speed",
    3: "Projectile mass",
    4: "Projectile cross-sectional area",
    5: "Projectile drag coefficient",
    6: "Planet mass",
    7: "Atmospheric fluid density",
    8: "Planet radius"}

    unitMap = {1: "/ s", 2: "/ m/s", 3: "/ kg", 4: "/ m^2", 5: "", 6: "/ kg", 7: "/ kg/m^3", 8: "/ m"}

    plt.xlabel("{} {}".format(stringMap.get(variable), unitMap.get(variable)))
    plt.ylabel("Optimal Angle / Degrees")
    plt.title("Relationship Between {} and Optimal Angle Under Given Conditions".format(stringMap[variable]))
    plt.grid(True)
    plt.plot(xValues,yValues)
    plt.show()

if __name__ == "__main__":
    main()

