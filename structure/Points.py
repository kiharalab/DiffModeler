import numpy as np
from structure.acc_utils import acc_merge_point
class Points(object):
    def __init__(self,params,cnt):
        self.params=params
        self.cd=np.zeros([cnt,3])#Coordinates record
        self.dens=np.zeros(cnt)#Density records
        self.origrid=np.zeros([cnt,3])#gridline record
        self.ori_dens=np.zeros(cnt)
        self.Ncd=0#Number of effective coordinates
        self.Nori=0#Number of effective grid points

    def recover_density(self):
        self.dens = self.dens*self.Nori

    def Merge_point(self, mrc, point_path):
        dcut = self.params['m'] / mrc.widthx#distance cut off
        d2cut = dcut ** 2
        rdcut = self.params['f']
        if self.Ncd == 0:
            print('not running merge shifting parts')
            return
        # id number of each points, so that when i is merged to j,
        # we can mark i with j, so that we know where it goes
        self.member = np.zeros(self.Ncd)
        self.stock = np.ones(self.Ncd)  # label for each points
        # Find the dmin and dmax in
        dmin = np.min(self.dens)
        dmax = np.max(self.dens)
        for i in range(self.Ncd):
            self.member[i] = i
        drange = dmax - dmin
        print('here we get the density range %f' % drange)
        rv_range = 1.0 / drange
        tmp_Ncd = int(self.Ncd)
        tmp_dens = np.array(self.dens)
        tmp_stock = np.array(self.stock)
        tmp_cd = np.array(self.cd)
        tmp_member = np.array(self.member)

        self.stock, self.member = acc_merge_point(tmp_Ncd, tmp_dens, dmin,
                                                  rv_range, rdcut, tmp_stock,
                                                  tmp_cd, d2cut,tmp_member)

        # Make new save arrays for
        Nmerge = 0
        # Get now the number of points
        for i in range(self.Ncd):
            if self.stock[i] == 1:
                Nmerge += 1

        self.merged_data = np.zeros([self.Ncd, 6])
        # init_id, x,y,z,density, merged_to_id
        self.merged_cd_dens = np.zeros([Nmerge, 4])  # first 3 column is the coordinates, the last 1 is the density
        # 1st column is the id points to the merged points, which is the member id
        cnt = 0  # special count for the left merged data
        self.Nmerge = Nmerge
        for i in range(self.Ncd):
            if self.stock[i]:
                self.merged_data[i][5] = cnt
                for k in range(3):
                    self.merged_cd_dens[cnt, k] = self.cd[i, k]
                self.merged_cd_dens[cnt][3] = self.dens[i]
                cnt += 1
            else:
                self.merged_data[i, 5] = -1  # indicates that point has been merged
            self.merged_data[i, 0] = self.member[i]
            for k in range(3):
                self.merged_data[i, k + 1] = self.cd[i, k]
            self.merged_data[i, 4] = self.dens[i]
        print('merging finishing with %d left' % self.Nmerge)
        np.savetxt(point_path, self.merged_data)
        np.savetxt(point_path[:-4] + 'onlymerged.txt', self.merged_cd_dens)
        self.normalize()

    def load_merge(self, mrc, point_path):
        print('load merge results')
        if self.Ncd == 0:
            print('not running merge shifting parts')
            return

        self.merged_data = np.loadtxt(point_path)
        self.merged_cd_dens = np.loadtxt(point_path[:-4] + 'onlymerged.txt')
        self.Nmerge = len(self.merged_cd_dens)

        print('merging finishing with %d left' % self.Nmerge)
        self.normalize()

    def normalize(self):
        self.merged_cd_dens[:, 3] = (self.merged_cd_dens[:, 3] - np.min(self.merged_data[:, 4])) / (
                    np.max(self.merged_data[:, 4]) - np.min(self.merged_data[:, 4]))+0.01#for 0 happen, which will raise problems in edge density calculation
        # self.min_dens=np.min(self.merged_data[:,4])
        # self.max_dens=np.max(self.merged_data[:,4])
        self.merged_data[:, 4] = (self.merged_data[:, 4] - np.min(self.merged_data[:, 4])) / (
                    np.max(self.merged_data[:, 4]) - np.min(self.merged_data[:, 4]))

    def clean_isolate_data(self,graph):
        self.mask = np.zeros(self.Nori)
        tmp_mask = np.zeros(self.Nori)
        for i in range(graph.Ne):
            if graph.edge[i].mst_label == False:
                continue
            v1 = graph.edge[i].id1
            v2 = graph.edge[i].id2
            tmp_mask[v1] = 1.00  # Member
            tmp_mask[v2] = 1.00  # Member of the tree
        for ii in range(self.Nori):
            m1 = int(self.merged_data[ii, 0])
            merged_id = int(self.merged_data[m1, 5])  # get the point id after merged
            if merged_id == -1:
                # print('not merged single point')
                continue
            self.mask[ii] = tmp_mask[merged_id]  # mark this original point that it is in edge or not
        print("We have %d/%d isolated points"%(len(np.argwhere(self.mask==0)),self.Nori))





