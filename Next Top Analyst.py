import sys
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
import numpy as np
import math
from itertools import product
from shapely.geometry import LineString, Point
from scipy import stats

#Plots a map showing where next top analyst could be found
#Rules dictate it is likely to be close to the river Spree, close to the Brandenburg Gate, and close to the satellite path

def viewProbHeatMaps(latAxisStart, latAxisStop, latAxisInterval, lonAxisStart, lonAxisStop, lonAxisInterval):
	#create grid of map points
	mapPoints = getAllPoints(latAxisStart, latAxisStop, latAxisInterval, lonAxisStart, lonAxisStop, lonAxisInterval)

	#for each grid point - get probability analyst is at this location
	#caculate probability according to river rule
	probFromRiver = riverRule(mapPoints, 0, 2730/2)

	#calculate probability according to Barangenburg Gate rule
	probFromGate = gateRule(mapPoints, 4700, 3877, 52.516288, 13.377689)

	#calculate probability according to satellite rule
	probFromSatellite = satelliteRule(mapPoints, 52.590117, 13.39915, 52.437385, 13.553989, 6371, 2400/2, 0)

	#plot heatmap of each probability
	#plotHeatMap(probFromRiver, mapPoints)
	plotHeatMap(probFromGate, mapPoints)
	#plotHeatMap(probFromSatellite, mapPoints)

	#compile all probabilities
	probTotal = compileProbs(probFromRiver, probFromGate)

	#plot heatmap of all probabilities
	#plotHeatMap(probTotal, mapPoints)

def getAllPoints(latStart, latStop, latInterval, lonStart, lonStop, lonInterval):
	mapPointsLon = np.linspace(lonStart, lonStop, num=lonInterval)
	mapPointsLat = np.linspace(latStart, latStop, num=latInterval)
	mapPointsLonAll, mapPointsLatAll = np.meshgrid(mapPointsLon, mapPointsLat) 
	allPoints = []
	for i in range(0, len(mapPointsLonAll)):
		for j in range(0, len(mapPointsLonAll[i])):
			allPoints.append([mapPointsLonAll[i][j],mapPointsLatAll[i][j]])
	return allPoints

def riverRule(allPoints, mu, sigma):
	#get river coordinates
	riverCoords = getRiverCoords()

	#loop through map to get probabilities of analyst being at location - due to river
	probRivers=[]
	for coords in allPoints:
		iLon, iLat = getLonLat(coords)

		#find shortest distance from river
		shortestDis = getShortestDisLinePoint(iLat, iLon, riverCoords)

		#for distance from river, get prob of analyst existence
		probRiver = mlab.normpdf(shortestDis,mu,sigma)
		probRivers.append(probRiver)
	return probRivers

def getRiverCoords():
	#parse txt file and create list of river GPS coords
	riverFile = open(sys.argv[1])
	riverCoords = []
	for line in riverFile:
		lat, lon = line.split(",")
		riverCoords.append([float(lon), float(lat)])
	return riverCoords

def getLonLat(coordPair):
	lon = coordPair[0]
	lat = coordPair[1]
	return lon, lat

#SARAH something is up with this distance, play with it
def getShortestDisLinePoint(pointLat, pointLon, lineCoords):
	#et shortest distance between point and line
	line = LineString(lineCoords) 
	point = Point(pointLon, pointLat)
	distance = point.distance(line)
	return distance

def gateRule(allPoints, mu, mode, fixedPointLat, fixedPointLon):
	#find sd using mode(X)=exp(mu-var)
	var = mu-math.log(mode)
	sigma = math.sqrt(var)

	muLog = np.log(mu)
	sigmaLog = np.log(sigma)
	
	#loop through map to get probabilities of analyst being at location - due to gate
	probGates=[]
	for coords in allPoints:
		iLon, iLat = getLonLat(coords)

		#find shortest distance from gate
		shortestDis = getShortestDisPointPoint(iLat, iLon, fixedPointLat, fixedPointLon)

		#for distance from gate, get prob of analyst existence
		if shortestDis > 0:
			probGate = stats.lognorm.pdf(shortestDis, sigmaLog, 0, np.exp(muLog))
		else:
			probGate = 0

		probGates.append(probGate)
	print stats.lognorm.pdf(0, sigmaLog, 0, np.exp(muLog))
	print stats.lognorm.pdf(0.001, sigmaLog, 0, np.exp(muLog))
	return probGates

def getShortestDisPointPoint(pointLat1, pointLon1, pointLat2, pointLon2): 
	point1 = Point(pointLat1, pointLon1) 
	point2 = Point(pointLat2, pointLon2)
	distance = point1.distance(point2)
	return distance

def lognstat(mu, var):
    """Calculate the mean of and variance of the lognormal distribution given
    the mean (`mu`) and standard deviation (`sigma`), of the associated normal 
    distribution."""
    m = np.exp(mu + (var / 2.0))
    v = np.exp(2 * (mu + var)) * (np.exp(var)-1)
    return m, v

def satelliteRule(allPoints, latStart, lonStart, latEnd, lonEnd, earthR, sigma, mu):
	#get satellite path coordinates
	pathCoords = ([lonStart, latStart], [lonEnd, latEnd])

	#loop through map to get probabilities of analyst being at location - due to river
	probPaths=[]
	for coords in allPoints:
		iLon, iLat = getLonLat(coords)

		#find shortest distance from river
		shortestDis = getShortestDisLinePoint(iLat, iLon, pathCoords)

		#for distance from river, get prob of analyst existence
		probPath = mlab.normpdf(shortestDis,mu,sigma)
		probPaths.append(probPath)
	return probPaths


def compileProbs(probs1, probs2):
	totalProbs = []
	for i in range(0, len(probs1)):
		totalProb = probs1[i] * probs2[i]
		totalProbs.append(totalProb)
	return totalProbs

def plotHeatMap(prob, allPoints):
	pyplot.scatter(list(zip(*allPoints)[0]), list(zip(*allPoints)[1]), c=prob, s=100, edgecolors='none', marker="h")
	pyplot.colorbar().set_label('probability of finding top analyst')
	pyplot.show()

def main():
	#overview with whole river spree
	viewProbHeatMaps(52.00, 53.00, 100, 12.80, 14.20, 100)

	#focus down to bridge area
	viewProbHeatMaps(52.20, 52.70, 100, 13.10, 13.60, 100)

if __name__ == '__main__':
	main()