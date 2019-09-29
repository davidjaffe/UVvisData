#!/usr/bin/env python
'''
simulation of the 1ton prototype
units mm, MeV, nm
20140521
'''
import math
import sys
#import random

import datetime
import numpy
import copy

class climate():
    def __init__(self):

        self.debug = 1

        self.GADfile = 'CDATA/GlobalAirportDatabase/GlobalAirportDatabase.txt'
        self.Airfarefile = 'CDATA/Airfares_CityPairs_20190929.csv'

        self.originalTeamList = ['Arizona Diamondbacks',
                          'Atlanta Braves',
                          'Baltimore Orioles',
                          'Boston Red Sox',
                          'Chicago Cubs',
                          'Chicago White Sox',
                          'Cincinnati Reds',
                          'Cleveland Indians',
                          'Colorado Rockies',
                          'Detroit Tigers',
                          'Miami Marlins',
                          'Houston Astros',
                          'Kansas City Royals',
                          'Los Angeles Angels of Anaheim',
                          'Los Angeles Dodgers',
                          'Milwaukee Brewers',
                          'Minnesota Twins',
                          'New York Mets',
                          'New York Yankees',
                          'Oakland Athletics',
                          'Philadelphia Phillies',
                          'Pittsburgh Pirates',
                          'St. Louis Cardinals',
                          'San Diego Padres',
                          'San Francisco Giants',
                          'Seattle Mariners',
                          'Tampa Bay Rays',
                          'Texas Rangers',
                          'Toronto Blue Jays',
                          'Washington Nationals']
        # team names altered to serve as keys and team airport city
        self.teams = {'Arizona_Diamondbacks':'Phoenix',
                          'Atlanta_Braves':'Atlanta',
                          'Baltimore_Orioles':'Baltimore',
                          'Boston_Red_Sox':'Boston',
                          'Chicago_Cubs':'Chicago',
                          'Chicago_White_Sox':'Chicago',
                          'Cincinnati_Reds':'Cincinnati',
                          'Cleveland_Indians':'Cleveland',
                          'Colorado_Rockies':'Denver',
                          'Detroit_Tigers':'Detroit',
                          'Miami_Marlins':'Miami',
                          'Houston_Astros':'Houston',
                          'Kansas_City_Royals':'Kansas City',
                          'Los_Angeles_Angels_of_Anaheim':'SANTA ANA',
                          'Los_Angeles_Dodgers':'Los Angeles',
                          'Milwaukee_Brewers':'Milwaukee',
                          'Minnesota_Twins':'Minneapolis',
                          'New_York_Mets':'New York',
                          'New_York_Yankees':'New York',
                          'Oakland_Athletics':'Oakland',
                          'Philadelphia_Phillies':'Philadelphia',
                          'Pittsburgh_Pirates':'PITTSBURGH (PENNSYLVA)',
                          'St._Louis_Cardinals':'St. Louis',
                          'San_Diego_Padres':'San Diego',
                          'San_Francisco_Giants':'San Francisco',
                          'Seattle_Mariners':'Seattle',
                          'Tampa_Bay_Rays':'Tampa',
                          'Texas_Rangers':'DALLAS-FORT WORTH',
                          'Toronto_Blue Jays':'Toronto',
                          'Washington_Nationals':'Washington'}
        
        # sort keys that are team names to be alphabetical
        self.teamList = self.teams.keys()
        self.teamList.sort()

        self.position = None
        
        
        return
    def readGlobalAirportDatabase(self):
        '''
        gets latitude, longitude of cities with MLB teams
        Data retrieved 20190915 from partow.net/miscellaneous/airportdatabase/#Download
        Field	Name	Type
01	ICAO Code	String (3-4 chars, A - Z)
02	IATA Code	String (3 chars, A - Z)
03	Airport Name	String
04	City/Town	String
05	Country	String
06	Latitude Degrees	Integer [0,360]
07	Latitude Minutes	Integer [0,60]
08	Latitude Seconds	Integer [0,60]
09	Latitude Direction	Char (N or S)
10	Longitude Degrees	Integer [0,360]
11	Longitude Minutes	Integer [0,60]
12	Longitude Seconds	Integer [0,60]
13	Longitude Direction	Char (E or W)
14	Altitude	Integer [-99999,+99999]
(Altitude in meters from mean sea level)
15	Latitude Decimal Degrees	Floating point [-90,90]
16	Longitude Decimal Degrees	Floating point [-180,180]
'''
        f = open(self.GADfile,'r')
        self.position = {}
        for line in f:
            s = line[:-1].split(':') # remove \n
            if self.debug>2: print 'climate.readGAD s',s
            code,airport,city,country,latitude,longitude = s[1],s[2],s[3],s[4],float(s[14]),float(s[15])
            if s[1]!='N/A': # Two Washington airports have lat,long=0,0
                names = self.mlbTeam(city,country)
                if names is not None:
                    for team in names:
                        self.position[team] = (latitude,longitude)
        f.close()
        print 'readGlobalAirportDatabase Found',len(self.position),'mlb cities'
        return
    def getGAD(self,reportDup=False):
        '''
        return dict of Global Airport Data in form dict[XXX] where XXX is three letter code IATA of airport

        GAD validity checks:
        IATA code not N/A
        city not N/A
        IATA code matches last 3 letters in IACO code (Akron and Akure, Nigeria have same IATA code AKR)

        set reportDup = True to report duplicate IATA
        '''
        f = open(self.GADfile,'r')
        GAD = {}
        for line in f:
            s = line[:-1].split(':')
            IACO = s[0]
            IATA = s[1]
            city = s[2]
            if IATA!='N/A' and city!='N/A':

                if IATA in GAD:   # already in dict, is this the best match?
                    if IATA==IACO[1:]: # this is better match
                        GAD[IATA] = s
                    elif IATA==GAD[IATA][0][1:]: # already have best match
                        pass
                    else:
                        if reportDup: print 'climate.readGAD DUPLICATE IATA',IATA,'in line',line[:-1]
                        #sys.exit('climate.readGAD ERROR IATA '+IATA+' IACO '+IACO)
                else:
                    GAD[IATA] = s
        print 'climate.getGAD Processed',self.GADfile
        f.close()
        return GAD
    def findIATA(self,GAD,city,country=None,debug=0):
        '''
        given dict of Global Airport Data
        return IATA code corresponding to input city, country
        else return None
        '''
        icity = 3
        icountry = 4
        ucity = city.upper()
        ucountry = None
        if country is not None: ucountry = country.upper()
        for IATA in GAD:
            CITY,COUNTRY = GAD[IATA][icity],GAD[IATA][icountry]
            if ucountry is None: COUNTRY = None
            if debug>0: print 'climate.findIATA ucity',ucity,'ucountry',ucountry,'CITY',CITY,'COUNTRY',COUNTRY
            if CITY==ucity and COUNTRY==ucountry : return IATA
        return None
    def readAirFares(self):
        '''
        Return dict[city1] = [fare,city1,city2] of airfares for city1,city2 pairs where city2==Tokyo 
        from airfare file
        0 = city1
        1 = city2
        3 = alternate name for city1 (eg., Bombay for Mumbai)
        2 = fare in USD
        4 = comments (airline nonstop or onestop)
        '''
        f = open(self.Airfarefile,'r')
        AirFares = {}
        for line in f:
            if 'Origin' not in line:
                s = line[:-1].split(',')
                city1 = s[0]
                city2 = s[1]
                acity1= s[3]
                fare  = float(s[2])
                if city1 in AirFares:
                    print 'climate.readAirFares ERROR line',line[:-1]
                    sys.exit('climate.readAirFares ERROR Duplicate city1 '+city1)
                else:
                    AirFares[city1] = [fare,city1,city2,acity1]
        f.close()
        print 'self.readAirFares Processed',self.Airfarefile
        return AirFares
                    
    def mlbTeam(self,city,country):
        '''
        return list of MLB team names if city,country combination is MLB city, otherwise None
        '''
        names = []
        if country=='USA' or country=='CANADA':
            for team in self.teams:
                if city.lower() == self.teams[team].lower():
                    if self.debug>1:
                        print 'climate.mlbTeam Match',city,country,'to',team,self.teams[team]
                    names.append(team)
        if len(names)>0: return names
        return None
    def getTeamPosition(self):
        '''
        assign longitude,latitude to each team
        '''
        self.readGlobalAirportDatabase()
        print '\nclimate.getTeamPosition: report latitude,longitude for airport near each team home'
        for team in self.teamList:
            if team in self.position:
                print team,self.position[team],
            else:
                print team,'**** NO POSITION *****'
        print '\n'
        return
    def getDistances(self):
        '''
        get distances between all pairs of teams
        return dict of distances of all pairs of teams
        '''
        pairDistance = {}
        n = 0
        if self.debug>1: print '\nDistance between all teams'
        
        for i1,team1 in enumerate(self.teamList):
            for team2 in self.teamList[i1+1:]:
                distance = self.haversine( self.position[team1],self.position[team2] )
                pairDistance[n] = [distance,team1,team2]
                if self.debug>1: print n,distance,team1,team2
                n += 1
        spD = sorted(pairDistance.items(), key=lambda x:x[1])
        print 'climate.getDistances Nearest',spD[0],'Farthest',spD[-1]
        return pairDistance
    def reDistribute(self):
        '''
        investigate schemes to redistribute teams to divisions
        '''
        pairDistance = self.getDistances() # pairDistance[n] = [distance,team1,team2]
        cities = ['Seattle','Miami','Chicago','Boston']
        anchors = []
        for city in cities:
            for team in self.teamList:
                if city in team:
                    anchors.append(team)
                    break
        print '\nclimate.reDistribute anchor teams',anchors
        for anchor in anchors:
            pD = {}
            for n in pairDistance:
                d,t1,t2 = pairDistance[n]
                if anchor==t1: pD[t2] = d
                if anchor==t2: pD[t1] = d
            spD = sorted(pD.items(), key=lambda x:x[1])
            print '\nclimate.reDistribute anchor',anchor,'\nteams',spD[:8]
        return
        
    def haversine(self,coord1, coord2):
        R = 6372800  # Earth radius in meters
        R = R/1000. # radius in km
        lat1, lon1 = coord1
        lat2, lon2 = coord2

        phi1, phi2 = math.radians(lat1), math.radians(lat2) 
        dphi       = math.radians(lat2 - lat1)
        dlambda    = math.radians(lon2 - lon1)

        a = math.sin(dphi/2)**2 + \
            math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2

        return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))
    def test(self):
        london_coord = 51.5073219,  -0.1276474
        cities = {
            'berlin': (52.5170365,  13.3888599),
            'vienna': (48.2083537,  16.3725042),
            'sydney': (-33.8548157, 151.2164539),
            'madrid': (40.4167047,  -3.7035825) 
            }

        for city, coord in cities.items():
            distance = self.haversine(london_coord, coord)
            print(city, distance)
        return
    def mainAirFares(self):
        '''
        main routine to figure out total of airfares for B2GM

        First get airport data, then airfares between city pairs (2d city is always Tokyo),
        then compute distances and fare/km
        '''
        GAD = self.getGAD()
        AirFares = self.readAirFares()
        cityIATA = {}
        Nfare,avcpkm,mincpkm,maxcpkm,rms = 0, 0., 1.e20, -1.e20, 0.
        for city1 in AirFares:
            fare = AirFares[city1][0]
            city2= AirFares[city1][2]
            acity1= AirFares[city1][3] # alternate name of city1
            debug,country = 0,None
            IATA = ct.findIATA(GAD,city1,country=country,debug=debug)
            if IATA is None: IATA = ct.findIATA(GAD,acity1,country=country,debug=debug)
            if IATA is None:
                sys.exit('climate.mainAirFares ERROR No IATA for '+city1)
            else:
                cityIATA[city1] = [IATA,float(GAD[IATA][14]),float(GAD[IATA][15])] # IATA, Latitude, Longitude
                
            if city2 not in cityIATA: # find city2 in database
                IATA = ct.findIATA(GAD,city2)
                if IATA is None:
                    sys.exit('climate.mainAirFares ERROR No IATA for '+city2)
                else:
                    cityIATA[city2] = [IATA,float(GAD[IATA][14]),float(GAD[IATA][15])]
            p1 = cityIATA[city1][1:]
            p2 = cityIATA[city2][1:]
            distance = self.haversine(p1,p2)
            cpkm = None
            if distance>0:cpkm = fare/distance
            Nfare += 1
            avcpkm += cpkm
            rms += cpkm*cpkm
            maxcpkm = max(maxcpkm,cpkm)
            mincpkm = min(mincpkm,cpkm)
            print 'self.mainAirFares',city1,city2,'Fare(USD)',fare,'distance(km)',distance,'USD/km',cpkm
        avcpkm = avcpkm/float(Nfare)
        rms = math.sqrt( float(Nfare)/float(Nfare-1) * (rms/float(Nfare) - avcpkm*avcpkm) )
        print 'self.mainAirFares #',Nfare,'average,sigma,min,max(USD/km)',avcpkm,rms,mincpkm,maxcpkm
        return
if __name__ == '__main__' :
   
    ct = climate()
    MLB = False
    if MLB:
        ct.getTeamPosition()
        ct.reDistribute()
    else:
        ct.mainAirFares()
            
#    ct.test()
