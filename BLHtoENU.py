from os import pardir
from pprint import pprint
import numpy as np
import math 
import csv
import pprint

#in this class, handle angle values as radians after inputting
class trans_coordinate():

    def __init__(self, data_name):
        self.name = data_name
        self.beta = 0
        self.lambda_ = 0
        self.ht = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.xyz = []
        self.ref = []
        self.enu = []
        
    # assume dd.dddd~ #transform to radians
    def input_BLH_deg(self, latitude, longitude, height):
        self.beta = math.radians(latitude)       #unit[deg] #beta #covert to radians
        self.lambda_ = math.radians(longitude)    #unit[deg] #lambda #covert to radians
        self.ht = height        #ellipsoidal height[m]


    # assume ddmmss.ssss~ #transform to radians
    def input_BLH_ddmmss(self, latitude, longitude, height):
        def ddmmss2deg(_x):
            dd = float((str(math.floor(_x)))[:-4])
            mm = float((str(math.floor(_x)))[-4:-2])
            ss = float((str(math.floor(_x)))[-2:]) + (_x - math.floor(_x))
            return dd + mm/60 + ss/60/60
        self.beta = math.radians(ddmmss2deg(latitude))     #beta #covert to radians
        self.lambda_ = math.radians(ddmmss2deg(longitude))    #lambda #covert to radians
        self.ht = height        #ellipsoidal height[m]


    def BLHtoECEF(self):
        f = 1/298.257223563     #the flattening factor of the Earth
        a = 6378137.0    #semimajor axis
        b = a*(1.0-f)       #semiminor axis
        e = math.sqrt(2*f-f*f)      #eccentricity

        n = a/math.sqrt(1-e*e*math.sin(self.beta)*math.sin(self.beta))  #prime vertical

        self.x = (n + self.ht)*math.cos(self.beta)*math.cos(self.lambda_)
        self.y = (n + self.ht)*math.cos(self.beta)*math.sin(self.lambda_)
        self.z = (n*(1 - e*e) + self.ht)*math.sin(self.beta)
        self.xyz = np.matrix([self.x, self.y, self.z])    #ECEF x,y,z output
        return self.xyz



    def rot_x(self, sita):
        s = np.sin(sita)
        c = np.cos(sita)
        R_x = np.matrix((
            (1, 0, 0),
            (0, c, s),
            (0, -s, c)
        ))
        return R_x

    def rot_y(self, sita):
        s = np.sin(sita)
        c = np.cos(sita)
        R_y = np.matrix((
            (c, 0, -s),
            (0, 1, 0),
            (s, 0, c)
        ))
        return R_y

    def rot_z(self, sita):
        s = np.sin(sita)
        c = np.cos(sita)
        R_z = np.matrix((
            (c, s, 0),
            (-s, c, 0),
            (0, 0, 1)
        ))
        return R_z
        

    def ECEFtoENU(self, ref_beta, ref_lambda_, ref_xyz):
        self.ref = np.matrix(ref_xyz)      #
        #print((self.xyz-self.ref).T)
        
        enu_np_mat = self.rot_z(math.pi/2)*self.rot_y(math.pi/2-ref_beta)*\
            self.rot_z(ref_lambda_)*((self.xyz-self.ref).T)
        
        print(enu_np_mat)
        
        #return enu_np_mat
        
   


if __name__ == '__main__':
    path = "position.txt"
    pos_blh_info = []
    tmp = []

    with open(path) as f:
        reader = csv.reader(f)
        row_index = 1
        for row in reader:
           
            if(row_index%2 == 1 ):
                tmp.append(row)
                               
            if(row_index%10 == 0 and row_index != 0):
                
                pos_blh_info.append(tmp)
                tmp = []

            row_index  += 1
    pprint.pprint(pos_blh_info)
    labels = []
   
    for i in range(len(pos_blh_info)):
        labels.append(("".join(pos_blh_info[i][0]))[:-4])
        #print(("".join(pos_blh_info[i][0]))[:-4])

    
    instances = [0]*len(labels)
    for i in range(len(pos_blh_info)):
        instances[i] = trans_coordinate(labels[i])
        lat = float(pos_blh_info[i][1][1])
        lon = float(pos_blh_info[i][2][1])
        ht = float(pos_blh_info[i][3][1])
        instances[i].input_BLH_deg(lat,lon,ht)
        instances[i].BLHtoECEF()
        print(instances[i].name)
        print(instances[i].xyz)
        
    x = 1   #choose ENU transformation reference
    ref = instances[x]
    del instances[x]

    for i in range(len(instances)):
        print(instances[i].name)
        instances[i].ECEFtoENU(ref.beta, ref.lambda_, ref.BLHtoECEF())

