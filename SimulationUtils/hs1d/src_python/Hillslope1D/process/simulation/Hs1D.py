import numpy as np
#from scipy.spatial.distance import cdist
#from scipy.interpolate import UnivariateSpline

class Hs1D(object):
    """
        Class containing the hillslope properties. It's an attribute of
        SpaceDiscretization Class
        #######################################################################
        @☺param
            xmin : minimal coordinate of the hillslope (in meters)
            xmax : maximal coordinate of the hillslope (in meters)
            nx : number of meshes used to describe the hillslope
            angle : slope of the hillslope on each mesh
            w : width of the hillslope on each mesh
            soil_depth : depth of the reservoir on each mesh
            k : hydraulic conductivity on the hillslope (in m/s)
            f : kinematic porosity on the hillslope

        @attributes:
          k : hydraulic conductivity
          f : porosity
          angle_edges : slope on edges
          soil_depth_edges : thickness of the layer on edges
          w_edges : width on edges
        WARNING ! Each property is defined over the edges of the mesh
        #######################################################################
    """

    def __init__(self, nx, angle, w, soil_depth, k, f, z_custom, x):
        #Extracting angle, w, x and soil_depth

        if isinstance(z_custom,int):
            if isinstance(angle, int) or isinstance(angle, float):
                if angle == 0:
                    self.angle_edges = np.arctan(0.05)*np.ones((nx+1, 1))
                elif angle > 0:
                    self.angle_edges = angle*np.ones((nx+1, 1))
                elif angle < 0:
                    print('Error, angle must be positive. Using absolute value')
                    self.angle_edges = np.abs(angle)*np.ones((nx+1, 1))
            else:
                self.angle_edges = np.reshape(angle, (len(angle), 1))
        else:
            diff_x = np.diff(np.reshape(x,(1,len(x))))
            diff_x = np.append(diff_x,diff_x[-1][-1])
            diff_x = np.reshape(diff_x,(len(diff_x),1))

            diff_z = np.diff(np.reshape(z_custom,(1,len(z_custom))))
            diff_z = np.append(diff_z,diff_z[-1][-1])
            diff_z = np.reshape(diff_z,(len(diff_z),1))
            self.angle_edges = np.arctan(diff_z/diff_x)

        if isinstance(w, int) or isinstance(w, float):
            if w == 0:
                self.w_edges = np.ones((nx+1, 1))*500
            elif w > 0:
                self.w_edges = np.ones((nx+1, 1))*w
            elif w < 0:
                print('Error, width must be positive. Using absolute value')
                self.w_edges = np.ones((nx+1, 1))*np.abs(w)
#            self.w = np.ones((nx+1,1))*100
        else:
            self.w_edges = np.reshape(w, (len(w), 1))

        if isinstance(soil_depth, int) or isinstance(soil_depth, float):
            if soil_depth == 0:
                self.soil_depth_edges = np.ones((1, nx+1))
            elif soil_depth > 0:
                self.soil_depth_edges = np.ones((1, nx+1))*soil_depth
            elif soil_depth < 0:
                print('Error, cells thickness must be a positive value. Using ', \
                      'absolute value')
                self.soil_depth_edges = np.ones((1, nx+1))*np.abs(soil_depth)
        else:
            self.soil_depth_edges = np.reshape(soil_depth, (len(soil_depth),1))

#        Set Hydraulic Parameters
        k = k                              #Convert k from m/h to m/s
        self.k = k
        self.f = f

#    def set_parameters(self,distance,z):
#         #Establishing the grid based on distances
#        if 'distance' in dir():
#            self.x_edges = np.arange(0,max(distance)+1,100)
#            self.x_edges = self.x_edges.reshape(len(self.x_edges),1)
#        else:
#            self.x_edges = np.arange(self.nx+1).reshape(self.nx+1,1)
#
#        Distance = cdist(distance[np.newaxis,:],self.x_edges[np.newaxis,:])
#        min_pos  = np.argmin(Distance,axis = 1)
#        if not 'z' in dir():
#            z = np.arange(len(self.x))
#        if min_pos:
#            self.z_edges = np.mean(np.bincount(min_pos,weight=z))
#            self.w_edges = np.mean(np.bincount(min_pos,weight=1))
#        else:
#            print('Cannot set parameters for hillslope ',str(self.Id),'\n')
#        self.soil_depth = np.ones(len(self.x))
#
#    def transform_to_constant_slope(self):
#        parameters = np.polyfit(self.x_edges,self.z_edges,1)
#        self.z_pred_edges = parameters[1]*self.x_edges + parameters[2]
#        self.R2_edges = np.sum((self.z_pred_edges - np.mean(self.z_edges))**2) / np.sum((self.z_edges - np.mean(self.z_edges))**2)
#        self.angle_edges = np.arctan((self.z_pred_edges[-1] - self.z_pred_edges[0])/(self.x_edges[-1]-self.x_edges[0]))
#        self.angle_edges = np.ones(np.shape(self.x_edges))
#
#        return self.z_pred_edges, self.R2_edges, self.angle_edges
#
#    def transform_to_spline_slope(self,SmoothingParam):
#        if len(self.x_edges)>2:
#            f3 = UnivariateSpline(self.x_edges,self.z_edges,s=SmoothingParam)
#            self.z_pred_edges = f3(self.x_edges)
#            x2 = self.x_edges+1
#            z2 = f3(x2)
#            x3 = self.x_edges-1
#            z3 = f3(x3)
#            self.R2_edges = np.sum((self.z_pred_edges - np.mean(self.z_edges))**2) / np.sum((self.z_edges - np.mean(self.z_edges))**2)
#            self.angle_edges = (z2-z3)/2
#            return self.z_pred_edges, self.R2_edges, self.angle_edges
#
#    def compute_elevation(self):
#        x_mean = (self.x_edges[1:] + self.x_edges[0:-2])/2
#        angle_mean = (self.angle_edges[1:] + self.angle_edges[0:-2])/2
#        dx = (self.x_edges[1:] - self.x_edges[0:-2])
#        self.elevation_edges = np.cumsum(dx * np.tan(angle_mean))
#        self.z_edges = np.interp(self.x_edges,x_mean,self.elevation_edges)
#        self.s_edges = np.cumsum(dx/np.cos(angle_mean))
#        return self
#
#
#    def get_x_edges(self):
#        return self.x_edges

    def get_w_edges(self):
        return self.w_edges

    def get_soil_depth_edges(self):
        return self.soil_depth_edges

    def get_angle_edges(self):
        return self.angle_edges

    def get_k(self):
        return self.k

    def get_f(self):
        return self.f







