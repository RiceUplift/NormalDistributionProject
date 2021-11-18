import numpy as np

def unnormalizedDist(x, mu, sigma):
    return(np.exp(-.5*((x-mu)/sigma)**2))

def riemannSum(mu, sigma, dx, x0, N):
    area = 0
    x = x0 #starting value of x
    for i in range(0, N):
        area += unnormalizedDist(x, mu, sigma)*dx
        x += dx
    return(area)

def normalDist(x, mu, sigma, C): #Calculates the value of the normalized function P(x;mu,sigma)
    return(C*np.exp(-.5*((x-mu)/sigma)**2))

def areaOneStdDev(mu, sigma, N, C):
    area = 0
    x0 = mu - sigma
    xf = mu + sigma
    dx = (xf - x0)/N
    x = x0 #starting value of x
    for i in range(0, N):
        area += normalDist(x, mu, sigma, C)*dx
        x += dx
    return(area)

def nStdDevArea(mu, sigma, N, C, n):
    area = 0
    x0 = mu - n*sigma
    xf = mu + n*sigma
    dx = (xf - x0)/N
    x = x0 #starting value of x
    for i in range(0, N):
        area += normalDist(x, mu, sigma, C)*dx
        x += dx
    return(area)

def errorFunction(mu, sigma, N, C, a):
    area = 0
    dx = 2*a/N #(a - (-a) = 2a)
    x = -a #starting value of x
    for i in range(0, N):
        area += normalDist(x, mu, sigma, C)*dx
        x += dx
    return(area)

#~~~~~~~~~MAIN CODE BELOW HERE~~~~~~~~~#

mean = 0 #mean
stdDev = 1 #standard deviation

N= 1000000 #number of rectangles used to approximate the area
x0 = -50 #starting x-value for the area
xf = 50 #end x-value for the area. You don't have to do it this way, but I think it makes it easier.
deltax = (xf - x0)/N #equal width for all triangles between x0 and xf
print( riemannSum(mean, stdDev, deltax, x0, N) ) #calculates the normalization constant. Only needed to help find C
C = 1/riemannSum(mean, stdDev, deltax, x0, N)

print( areaOneStdDev(mean, stdDev, N, C)) #calculates the area one standard deviation from the mean

n=5
print( nStdDevArea(mean, stdDev, N, C, n) ) #calculates the area n standard deviations away. 
print( errorFunction(mean, stdDev, N, C, 2)) #calculates the area on (-a,a). As written, this calculates the area on (-2,2). 2 can be replaced with whatever a-value you want.

