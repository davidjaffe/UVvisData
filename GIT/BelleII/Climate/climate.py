#!/usr/bin/env python
'''
climate and airfare related calculations
201909xx
'''
import math
import sys
#import random

import datetime
import numpy
import copy

import matplotlib.pyplot as plt

class climate():
    def __init__(self):

        self.debug = 1

        self.GADfile = 'CDATA/GlobalAirportDatabase/GlobalAirportDatabase.txt'
        self.Airfarefile = 'CDATA/Airfares_CityPairs_20190929.csv'
        self.figdir = 'FIGURES/'

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
    def readB2I(self):
        '''
        read Belle II institutions lists to get inst names, country, city
        return dict[short name] = [city,country,long name]
        '''
        self.B2InstitutionsFile = 'CDATA/BelleII_institutions_20190929.csv'
        f = open(self.B2InstitutionsFile,'r')
        B2Inst = {}
        for line in f:
            if 'Postal' not in line: # avoids header
                s = line[:-1].split(',')
                sname = s[0]
                lname = s[1]
                country = s[4]
                city = s[5]
                B2Inst[sname] = [city,country,lname]
        f.close()
        print 'climate.readB2I Processed',len(B2Inst),'institutions in',self.B2InstitutionsFile
        return B2Inst
    def readB2M(self):
        '''
        read Belle II members to get Belle II ID, name and institution name
        return dict[B2id,firstname,middlename,lastname,lastnameprefix,inst]
        '''
        self.B2MembersFile = 'CDATA/BelleII_members_20190929.csv'
        f = open(self.B2MembersFile,'r')
        B2Members = {}
        for line in f:
            if 'B2id' not in line: # avoids header
                s = line[:-1].split(',')
                for i,x in enumerate(s):
                    y = x.replace('"','')
                    if y=='-': y = ''
                    s[i] = y
                B2id = s[0]
                firstn = s[1]
                lastn  = s[2]
                middlen= s[3]
                lastnpre=s[4]
                sname  = s[11] # institution short name
                B2Members[B2id] = [firstn,middlen,lastn,lastnpre,sname]
        f.close()
        print 'climate.ReadB2M Processed',len(B2Members),'members in',self.B2MembersFile
        return B2Members
    def readB2GM(self):
        '''
        return list of name1,name2 of cut/paste of B2GM participants https://indico.belle2.org/event/971/registrations/participants 
        as of 20190929
        Try to take into account that name1 or name2 can be the family name or given name
        '''
        self.B2GMfile = 'CDATA/B2GM_Oct2019_20190929.txt'
        debug = False
        f = open(self.B2GMfile,'r')
        B2GMfolks = []
        for line in f:
            if 'Participant' not in line: # avoid header
                j = line.find('\t')
                name1 = line[:j].replace('\t',' ').strip()
                name2 = line[j+1:-1].replace('\t',' ').strip()
                B2GMfolks.append( [name1,name2] )
        f.close()
        print 'climate.readB2GM Processed',len(B2GMfolks),'participants in',self.B2GMfile
        if debug: print B2GMfolks
        return B2GMfolks
    def drawIt(self,x,y,xtitle,ytitle,title,figDir=None,ylog=True,xlims=[200.,800.],ylims=[1.e-5,1.e-1]):
        '''
        draw graph defined by x,y

        '''
        debug = False
        if debug: print 'climate.drawIt len(x),len(y)',len(x),len(y),xtitle,ytitle,title,'xlimits',xlims,'ylimits',ylims
        if debug: print 'climate.drawIt x',x,'\ny',y
        plt.clf()
        plt.grid()
        plt.title(title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'
        if debug: print 'climate.drawIt figpdf',figpdf

        X = numpy.array(x)
        Y = numpy.array(y)
        if debug: print 'climate.drawit X',X,'\nY',Y
        plt.plot(X,Y,'o')
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : plt.yscale('log')
        if debug: print 'climate.drawit ylog',ylog
        plt.xlim(xlims)
        plt.ylim(ylims)

        if debug: print 'climate.drawit figDir',figDir
        
        if figDir is not None:
            figpdf = figDir + figpdf
            plt.savefig(figpdf)
            print 'climate.drawIt wrote',figpdf
        else:
            if debug: print 'climate.drawit show',title
            plt.show()
        return
    def goodMatch(self,a,b):
        return a.strip().lower()==b.strip().lower()
    def matchMtoI(self,Insts,Members,Parps):
        '''
        try to match participants Parps to Members and determine their institutions
        Result should be list of participants with their institutions including city, country
        '''
        debug = 1
        TooFew,TooMany,JustRight = 0,0,0
        for P in Parps:
            name1,name2 = P
            if debug>1: print 'Search for name1',name1,'name2',name2
            Matches = [] # list of B2id of potential matches
            for B2id in Members:
                firstn,middlen,lastn,lastnpre,sname = Members[B2id]
                if debug>1: print 'first,middle,last,pre',firstn,middlen,lastn,lastnpre
                if self.goodMatch(name1,lastn):
                    Matches.append(B2id)
                if self.goodMatch(name2,lastn):
                    if B2id not in Matches: Matches.append(B2id)

            if len(Matches)>1: # see if we can get a name1=firstn and name2=lastn or vice-versa
                newMatch = []
                for B2id in Matches:
                    firstn,middlen,lastn,lastnpre,sname = Members[B2id]
                    if self.goodMatch(name1,firstn) and self.goodMatch(name2,lastn): newMatch.append(B2id)
                    if self.goodMatch(name2,firstn) and self.goodMatch(name1,lastn) and B2id not in newMatch: newMatch.append(B2id)
                if len(newMatch)>0: Matches = newMatch
            L = len(Matches)
            
            if L==0:
                if debug>0: print name1,name2,'*** NO MATCHES *** '
                TooFew += 1
            elif L==1:
                JustRight += 1
            else:
                if debug>0:
                    print name1,name2,L,'matches',Matches
                    self.printMatches(Matches,Members)
                TooMany += 1
        print 'climate.matchMtoI',len(Parps),'participants',JustRight,'single matches',TooMany,'multi-matches',TooFew,'No matches'
        return
    def printMatches(self,idList,Members):
        for i,B2id in enumerate(idList):
            firstn,middlen,lastn,lastnpre,sname = Members[B2id]
            print ' match#',i,firstn,middlen,lastn,lastnpre
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
        distances,fares = [],[]
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
            fares.append(fare)
            distances.append(distance)
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
        print 'climate.mainAirFares #',Nfare,'average,sigma,min,max(USD/km)',avcpkm,rms,mincpkm,maxcpkm

        B2Inst = self.readB2I()
        B2Members = self.readB2M()
        B2GMfolks = self.readB2GM()
        self.matchMtoI(B2Inst,B2Members,B2GMfolks)
        
        
        if 0: 
            self.drawIt(distances,fares,'Distance between city and Tokyo in km','Fare(USD)','Cost per km',figDir=None,ylog=False,xlims=[0.,14000.],ylims=[0.,2400.])
        
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
