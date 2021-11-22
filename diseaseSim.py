# 
# diseaseSim.py - simulates the spread of disease through a population in a specified neighbourhood. Returns with the results of number of uninfected, infected, dead and immune people and plots a graph for its progress after each timesteps.
#
# 
#
# Version history:
#
# 25/4/19 - beta version released for FOP assignment
# 24/5/19 - extended version of the code
#
#

import numpy as np
import matplotlib.pyplot as plt
import random
import sys

if (len(sys.argv) != 13):                            #Runs default value if command line argument is less than 9 including the .py file
    print('\nRequires 12 variables for Population/Infected/Doctors/Number of Steps/p(infection)/p(death)/p(recovery)/lethal rate/Neighbourhood (M/V)/Row/Cols\n')
    print('\nUsing default values for initial population 1000, 15 infected, 20 Doctors, 80 steps, 0.5/0.001/0.4/0.05/0.1 probability of infected/death/recovery/immunity, 0.1 for lethal rate and Moore Neighbourhood\n')
    print('Dimension: 30 rows, 40 columns')
    INIT_POP = 1000
    INIT_INFECTED =15
    NUM_DOCS = 20
    NUM_STEPS = 80
    ProbInfected = 0.5
    ProbKilled = 0.001
    ProbRecovered = 0.25
    ProbImmune=0.05
    lethalrate = 0.1
    neighbchoice = 'M'
    NUM_ROWS=30
    NUM_COLS=40
else:                                              #Or else it will run the parameters that are typed in
    INIT_POP = int(sys.argv[1])
    INIT_INFECTED = int(sys.argv[2])
    NUM_DOCS = int(sys.argv[3])
    NUM_STEPS = int(sys.argv[4])
    ProbInfected = float(sys.argv[5])
    ProbKilled = float(sys.argv[6])
    ProbRecovered = float(sys.argv[7])
    ProbImmune = float(sys.argv[8])
    lethalrate= float(sys.argv[9])
    neighbchoice = str(sys.argv[10])
    NUM_ROWS = int(sys.argv[11])
    NUM_COLS = int(sys.argv[12])


def distuninf(grid, num_r, num_c, numpeep):      #Adds a random point in a matrix by '1' unless it's on the barriers for a given number of times
    for i in range(numpeep): 
        rpos = random.randint(1, (num_r-2))          
        while rpos == num_r//2:
              rpos = random.randint(1, (num_r-2))
        cpos = random.randint(1, (num_c-2))
        while cpos == num_c//2:
              cpos = random.randint(1, (num_c-2))
        grid[rpos, cpos] += 1 

def distinf(grid, num_r, num_c, numpeep):      #Adds a random point in a matrix by '1' at the left side of the grid also known as the infected area for a given number of times
    for i in range(numpeep):
        rpos = random.randint(1, (num_r-2))
        while rpos == num_r//2:
              rpos = random.randint(1, (num_r-2))
        cpos = random.randint(1, ((num_c-2)//2) - 1)
        grid[rpos, cpos] += 1

def distdoc(grid, num_r, num_c, numpeep):      #Adds a random point in a matrix by '1' in the bottom middle of the grid, this matrix represents the doctors dispatched
    for i in range(numpeep//2):
        rpos = random.randint(1, (num_r//2) - (num_r//3))
        cpos = random.randint(1, (num_c//2) - (num_c//3))
        grid[rpos, cpos] += 1
    for i in range(numpeep//2):
        rpos = random.randint(3*num_r//4, num_r-2 )
        cpos = random.randint(1 , (num_c//2) - (num_c//3))
        grid[rpos, cpos] += 1



def makeScatter(grid, num_r, num_c):                           #Makes a matrix suitable to plot on a scatterplot
    r_values = []
    c_values = []
    count_values = []
    for row in range(num_r):
        for col in range(num_c):
            if grid[row,col] > 0:
                r_values.append(row)
                c_values.append(col)
                count_values.append(grid[row,col]*20)
#                print("Value at (", row, ",", col, ") is ", grid[row, col])
    return(r_values, c_values, count_values)
    
def displayGrid(grid, num_r, num_c):                           #Displays the grid of a matrix
    for row in range(num_r-1, -1, -1):
        for col in range(num_c):
            print(grid[row,col], end=" ")
        print()

def barrier(grid, num_r, num_c):                               #Makes a barrier of '1' on the edges and in the middle of the column and row
    for row in range(num_r):                                 
         grid[row,0]= 1
         grid[row, num_c-1]= 1
         grid[row, num_c//2]= 1
    for col in range(num_c):
         grid[0, col]= 1
         grid[num_r-1, col]= 1
         grid[num_r//2, col]= 1
    grid[num_r-2, num_c//2]=0
    grid[1, num_c//2]=0
    grid[num_r//2 , 1]=0

def plotGrids():                                               #plot a scatterplot from the simulated data
    Irows, Icols, Icount = makeScatter(infected, NUM_ROWS, NUM_COLS)
    plt.scatter(Icols, Irows, s=Icount, c="r", alpha=0.5, label = "Infected")
    Urows, Ucols, Ucount = makeScatter(uninfected, NUM_ROWS, NUM_COLS)
    plt.scatter(Ucols, Urows, s=Ucount, c="b", alpha=0.5, label = "Uninfected")
    Drows, Dcols, Dcount = makeScatter(dead, NUM_ROWS, NUM_COLS)
    plt.scatter(Dcols, Drows, s=Dcount, c="k", alpha=0.5,label = "Dead")
    Hrows, Hcols, Hcount = makeScatter(doctor, NUM_ROWS, NUM_COLS)
    plt.scatter(Hcols, Hrows, s=Hcount, c="orange", alpha =1, label = "Doctor")
    Mrows, Mcols, Mcount = makeScatter(immune, NUM_ROWS, NUM_COLS)
    plt.scatter(Mcols, Mrows, s=Mcount, c="g", alpha=0.5, label = "Immune")
    Brows, Bcols, Bcount = makeScatter(world, NUM_ROWS, NUM_COLS)
    plt.scatter(Bcols, Brows, s=10, c ="k", marker = "x", alpha=1)
    plt.legend(prop= {'size':6}, loc = "upper right")
    plt.title("Disease Simulation no. "+str(timestep+1)+" Population: " +str(INIT_POP)+ " Infected: " +str(NumVar(infected))+" Doctors:"+str(NUM_DOCS))           
    plt.show()

def plotresults():                                                      #plots a linear graph that depicts the progress at the end of the simulation
    plt.plot(timestepcount,recordinfected, 'r-',label = "infected")
    plt.plot(timestepcount, recorduninfected, 'b-',label = "uninfected")
    plt.plot(timestepcount, recordimmune, 'g-', label = "immune")
    plt.plot(timestepcount, recorddead, 'k-', label = "dead")
    plt.title("Number of people after simulation. Pop/Inf/Doc: "+str(INIT_POP)+"/"+str(INIT_INFECTED)+"/"+str(NUM_DOCS)+"_"+str(neighbchoice))
    plt.ylabel('Number of people')
    plt.xlabel('Timesteps')
    plt.legend(prop= {'size':6}, loc = "upper right")
    plt.savefig("DiseaseSim_Pop:" + str(INIT_POP) + "_Inf"+ str(INIT_INFECTED) +"_Doctors:"+ str(NUM_DOCS)+".png")
    plt.show()

       
def movePeepsMooreStyle(cur, next, r, c):                      #Moves a person from one timestep to another in all 8 directions
#    print("Pos (", r, ",", c, ") has ", cur[r,c], " people")  
    for peep in range(cur[r,c]):                               
         rMove = random.randint(-1,1)
         cMove = random.randint(-1,1)
         if r == (NUM_ROWS-2) and c == (NUM_COLS//2):
            rMove = 0
         if r == 1 and c == (NUM_COLS//2):
            rMove = 0
         if r == ((NUM_ROWS)//2) and c == 1:
            cMove = 0
#         print("Move from (", r, ",", c, ") to (", r+rMove, "," , c+cMove, " 
         if (r + rMove) > (NUM_ROWS-2) or (r + rMove) < 1 or (r + rMove) == ((NUM_ROWS)//2):
             if (c + cMove) == 1:
                  if (r + rMove) == ((NUM_ROWS)//2):
                     rMove = rMove
                  else:
                     rMove = 0
             else:
                rMove = 0 
         if (c + cMove) > (NUM_COLS-2) or (c + cMove) < 1 or (c + cMove) == ((NUM_COLS//2)):
             if (r + rMove) == (NUM_ROWS -2) or (r + rMove) == 1:
                  if (c + cMove) == ((NUM_COLS//2)):
                     cMove = cMove
                  else:
                     cMove = 0
             else:
                cMove = 0   
         next[r + rMove, c + cMove] +=1
#         print("          (", r, ",", c, ") to (", r+rMove, "," , c+cMove, ")")


    
def movePeepsVonStyle(cur, next, r, c):                        #Moves a persom in up ,down ,left or right direction from one timestep to the next
    for peep in range(cur[r,c]):                               
         rMove = random.randint(-1,1)
         cMove = random.randint(-1,1)
         if r == (NUM_ROWS-2) and c == (NUM_COLS//2):
            rMove = 0
         if r == 1 and c ==(NUM_COLS//2):
            rMove = 0
         if r == ((NUM_ROWS)//2) and c == 1:
            cMove = 0
         if (rMove != 0):
             cMove = 0
         if (r + rMove) > (NUM_ROWS-2) or (r + rMove) < 1 or (r + rMove) == ((NUM_ROWS)//2):
             if (c + cMove) == 1:
                  if (r + rMove) == ((NUM_ROWS)//2):
                      rMove = rMove
                  else:
                      rMove = 0
             else:
                rMove = 0
         if (c + cMove) > (NUM_COLS-2) or (c + cMove) < 1 or (c + cMove) == ((NUM_COLS)//2):
             if (r + rMove) == (NUM_ROWS -2) or (r + rMove) == 1:
                  if (c + cMove) == ((NUM_COLS//2)):
                     cMove = cMove
                  else:
                     cMove = 0
             else:
                cMove = 0
         next[r + rMove,c + cMove] +=1
         #print("          (", r, ",", c, ") to (", r+rMove, "," , c+cMove, ")")

def staystill(cur, next, r, c):                                 #The person stays still when transitioning from one timestep to the next
    for peep in range(cur[r,c]):
         next[r ,c] +=1

def infect(inf, notinf, r, c, prob):                            #infect a non-infected person with a given probability
#    print("Pos (", r, ",", c, ") has ", inf[r,c], " inf people and ", notinf[r,c], " well people")
    prob = prob*inf[r,c]
    if prob:
        for peep in range(notinf[r,c]):
            if random.random() < prob:
                inf[r, c] +=1
                notinf[r, c] -=1
                #print("***** New infection (", r, ",", c, ")")

def kill(dead, inf, r, c, prob):                                 #Kills an infected person with a given probability
    prob = prob
    if prob:
        for peep in range(inf[r,c]):
            if random.random() < prob:
                dead[r, c] +=1
                inf[r,c] -=1
                #print("***** New Death (", r, ",",c, ")")

def recover(notinf, inf, immu, doc, r, c, prob, prob2):               #An infected person has a chance to recover.
    prob = prob + prob*immu[r,c] + doc[r,c]
    prob2 = prob2 + prob2*doc[r,c]
    if prob:
        for peep in range(inf[r,c]):
            if random.random() < prob:
                if random.random() < prob2:                      #If it recovers, it has a chance of becoming immune to the disease        
                    immu[r,c] +=1
                    inf[r,c] -=1
                    #print("***** Immune(", r, "," ,c, ")")
                else:
                    notinf[r,c] +=1
                    inf[r,c] -=1
                    #print("***** Recovered(",r,",",c,")")

def NumVar(variable):                                             #Converts an array matrix into a number of people in the selected category
    NonZeroVariables = variable[ np.nonzero(variable)]
    return sum(NonZeroVariables)

def recorddata():                                                  #Records data for each category after each for loop
     timestepcount[timestep+1] = timestep + 1
     recordinfected[timestep+1] = NumVar(infected)
     recorduninfected[timestep+1] = NumVar(uninfected)
     recordimmune[timestep+1] = NumVar(immune)
     recorddead[timestep+1] = NumVar(dead)


def shownumofppl():                                                #Displays the number of people in each category
    print("No. of Uninfected:", NumVar(uninfected))
    print("No. of Infected:", NumVar(infected))
    print("No. of Immune:", NumVar(immune))
    print("No. of People that died:", NumVar(dead))

def showstats():                                                #Shows results and probability of each
    print("\nInitially, there are", INIT_POP, "people, ", INIT_INFECTED, "are infected and", NUM_DOCS, " doctors are dispatched")
    print("Probability of Infection:", ProbInfected)
    print("Probability of Death:", ProbKilled)
    print("Probability of Recovery:",ProbRecovered)
    print("Probability of Immunity:", ProbImmune)
    print("Lethal Death Rate:", lethalrate)
    print("\nHere are the results after", NUM_STEPS, "timestep(s):\n")
    print("Death rate:", int((NumVar(dead)/INIT_POP)*100), "%")
    print("Immunity rate:", int((NumVar(immune)/INIT_POP)*100), "%")
    print("Infection rate:", int((NumVar(infected)/INIT_POP)*100), "%")

def showconclusion():                                                           #Analyses the death rate and infected rate and provides a conclusion
    deathrate = NumVar(dead)/INIT_POP
    infectedrate = NumVar(infected)/INIT_POP
    if deathrate > lethalrate :
       print("\nConclusion: EMERGENCY! ACTION REQUIRED IMMEDIATELY\n")
    else:
       if infectedrate > 0.15: 
           print("\nConclusion: Still too many infected in the population. Need more doctors\n")
       else:
          if NumVar(infected) == 0:
             print("\nConclusion: No more infected!\n")
          else:
             print("\nConclusion: Disease outbreak controlled. Enough doctors dispatched\n")
             
INIT_UNINFECTED = INIT_POP - INIT_INFECTED

world = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)                              #Generates zero matrices for each variable
infected = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
doctor = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
uninfected = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
dead = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
immune = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)

timestepcount = np.zeros(NUM_STEPS+1)                                             #Generates zero matrices to record the data for each time step
timestepcount[0] = 0
recordinfected = np.zeros(NUM_STEPS+1)
recordinfected[0] = INIT_INFECTED
recorduninfected = np.zeros(NUM_STEPS+1)
recorduninfected[0] = INIT_UNINFECTED
recordimmune = np.zeros(NUM_STEPS+1)
recorddead = np.zeros(NUM_STEPS+1)

barrier(world, NUM_ROWS, NUM_COLS)
distinf(infected, NUM_ROWS, NUM_COLS, INIT_INFECTED)
distuninf(uninfected, NUM_ROWS, NUM_COLS, INIT_UNINFECTED)
distdoc(doctor, NUM_ROWS, NUM_COLS, NUM_DOCS)
#print(world)
#print()
#displayGrid(infected, NUM_ROWS, NUM_COLS)
#print()
#displayGrid(uninfected, NUM_ROWS, NUM_COLS)
#displayGrid(world, NUM_ROWS, NUM_COLS)
timestep = -1
plotGrids()

neighbchoice = neighbchoice.upper()    
if neighbchoice == 'M':
   for timestep in range(NUM_STEPS):
       print("\n###################### TIMESTEP", timestep, "#####################\n")
       infected2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       uninfected2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       dead2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       immune2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       doctor2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       for row in range(NUM_ROWS):
           for col in range(NUM_COLS):
               infect(infected, uninfected, row, col, ProbInfected)
               kill(dead, infected, row, col, ProbKilled)
               recover(uninfected, infected, immune, doctor, row, col, ProbRecovered, ProbImmune)
               movePeepsMooreStyle(doctor, doctor2, row, col)
               movePeepsMooreStyle(infected, infected2, row, col)
               movePeepsMooreStyle(uninfected, uninfected2, row, col)
               staystill(dead, dead2, row, col)
               movePeepsMooreStyle(immune, immune2, row, col)
       recorddata()
       infected = infected2
       uninfected = uninfected2
       dead = dead2
       immune = immune2
       doctor = doctor2
       shownumofppl()
       plotGrids()
   print("\n####################Moore Neighbourhood selected######################\n")

if neighbchoice == 'V':
   for timestep in range(NUM_STEPS):
       print("\n###################### TIMESTEP", timestep, "#####################\n")
       infected2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       uninfected2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       dead2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       immune2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       doctor2 = np.zeros((NUM_ROWS, NUM_COLS), dtype=np.int)
       for row in range(NUM_ROWS):
           for col in range(NUM_COLS):
               infect(infected, uninfected, row, col, ProbInfected)
               kill(dead, infected, row, col, ProbKilled)
               recover(uninfected, infected, immune, doctor, row, col, ProbRecovered, ProbImmune)
               movePeepsVonStyle(infected, infected2, row, col)
               movePeepsVonStyle(uninfected, uninfected2, row, col)
               movePeepsVonStyle(immune, immune2, row, col)
               movePeepsVonStyle(doctor, doctor2, row, col)
               staystill(dead, dead2, row, col)
       recorddata()
       infected = infected2
       doctor = doctor2
       uninfected = uninfected2
       dead = dead2
       immune = immune2
       shownumofppl()
       plotGrids()
   print("\n####################Vonn Neumann Neighbourhood selected#################\n")

showstats()
shownumofppl()
showconclusion()
plotresults()
