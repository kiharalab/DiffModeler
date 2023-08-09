from numba import jit
import numpy as np
import math
@jit(nogil=True,nopython=True)
def carry_shift(point_cd,cnt,fmaxd,fsiv,xdim,ydim,zdim,dens):
    if True:
        point_dens=np.zeros(cnt)
        for i in range(cnt):
            if i%1000==0:
                print("mean shift",i,"/",cnt)
            #print(i)
            #print('start shifting for '+str(i))
            pos=np.zeros(3)
            for j in range(3):
                pos[j]=point_cd[i][j]
            if True:
                while True:
                    stp=np.zeros(3)
                    endp=np.zeros(3)
                    for j in range(3):
                        stp[j]=int(pos[j]-fmaxd)
                        if stp[j]<0:
                            stp[j]=0
                        endp[j]=int(pos[j]+fmaxd+1)
                    if endp[0]>=xdim:
                        endp[0]=xdim
                    if endp[1]>=ydim:
                        endp[1]=ydim
                    if endp[2]>=zdim:
                        endp[2]=zdim
                    dtotal=0
                    pos2=np.zeros(3)
                    for xp in range(int(stp[0]),int(endp[0])):
                        rx=float(xp-pos[0])**2
                        for yp in range(int(stp[1]),int(endp[1])):
                            ry=float(yp-pos[1])**2
                            for zp in range(int(stp[2]),int(endp[2])):
                                rz=float(zp-pos[2])**2
                                d2=rx+ry+rz#relative distance square
                                v=np.exp(-1.5*d2*fsiv)*dens[xp,yp,zp]#This is the bottom part of the equation, where pos represents y, (xp,yp,zp) represents xi
                                dtotal+=v
                                if v>0:
                                    pos2[0]+=v*(float)(xp)#pos2 is for the top part of the equation
                                    pos2[1]+=v*(float)(yp)
                                    pos2[2]+=v*(float)(zp)
                    if dtotal==0:
                        break
                    rd=1.00/float(dtotal)
                    tempcd=np.zeros(3)
                    for j in range(3):
                        pos2[j]*=rd#Now we get the equation result
                        tempcd[j]=pos[j]-pos2[j]
                        pos[j]=pos2[j]#Prepare for iteration
                    check_d=tempcd[0]**2+tempcd[1]**2+tempcd[2]**2#Iteration until you find the place is stable
                    if check_d<0.001:
                        break

            for j in range(3):

                point_cd[i][j]=pos[j]
            point_dens[i]=dtotal#/cnt
        return point_cd,point_dens

@jit(nogil=True,nopython=True)
def carry_shift_limit(point_cd,cnt,fmaxd,fsiv,xdim,ydim,zdim,dens,move_limit=3):
    if True:
        point_dens=np.zeros(cnt)
        for i in range(cnt):
            if i%1000==0:
                print(i)
            #print(i)
            #print('start shifting for '+str(i))
            pos=np.zeros(3)
            for j in range(3):
                pos[j]=point_cd[i][j]
            if True:
                while True:
                    stp=np.zeros(3)
                    endp=np.zeros(3)
                    for j in range(3):
                        stp[j]=int(pos[j]-fmaxd)
                        if stp[j]<0:
                            stp[j]=0
                        endp[j]=int(pos[j]+fmaxd+1)
                    if endp[0]>=xdim:
                        endp[0]=xdim
                    if endp[1]>=ydim:
                        endp[1]=ydim
                    if endp[2]>=zdim:
                        endp[2]=zdim
                    dtotal=0
                    pos2=np.zeros(3)
                    for xp in range(int(stp[0]),int(endp[0])):
                        rx=float(xp-pos[0])**2
                        for yp in range(int(stp[1]),int(endp[1])):
                            ry=float(yp-pos[1])**2
                            for zp in range(int(stp[2]),int(endp[2])):
                                rz=float(zp-pos[2])**2
                                d2=rx+ry+rz#relative distance square
                                v=np.exp(-1.5*d2*fsiv)*dens[xp,yp,zp]#This is the bottom part of the equation, where pos represents y, (xp,yp,zp) represents xi
                                dtotal+=v
                                if v>0:
                                    pos2[0]+=v*(float)(xp)#pos2 is for the top part of the equation
                                    pos2[1]+=v*(float)(yp)
                                    pos2[2]+=v*(float)(zp)
                    if dtotal==0:
                        break
                    rd=1.00/float(dtotal)
                    tempcd=np.zeros(3)
                    current_move = 0
                    for j in range(3):
                        pos2[j]*=rd#Now we get the equation result
                        tempcd[j]=pos[j]-pos2[j]
                        current_move += (pos2[j]-point_cd[i][j])**2
                        pos[j]=pos2[j]#Prepare for iteration
                    check_d=tempcd[0]**2+tempcd[1]**2+tempcd[2]**2#Iteration until you find the place is stable
                    if current_move>=move_limit**2:#limit further moving to avoid missing base predictions in some regions
                        break
                    if check_d<0.001:
                        break

            for j in range(3):

                point_cd[i][j]=pos[j]
            point_dens[i]=dtotal/cnt
        return point_cd,point_dens

@jit(nogil=True,nopython=True)
def acc_merge_point(Ncd,dens,dmin,rv_range,rdcut,stock,cd,d2cut,member):
    if True:
        for i in range(Ncd-1):
            if i%10000==0:
                print(i)
            tmp=np.zeros(3)
            if (dens[i]-dmin)*rv_range < rdcut:
                stock[i]=0#Label the small density parts as unused parts
            if stock[i]==0:
                continue
            for j in range(i+1,Ncd):
                if stock[j]==0:
                    continue
                d2=0
                for k in range(3):
                    tmp[k]=cd[i][k]-cd[j][k]
                    d2+=tmp[k]**2
                if d2<d2cut:
                    #Mark the merged points to where it goes
                    if dens[i]>dens[j]:
                        stock[j]=0
                        member[j]=i
                    else:
                        stock[i]=0
                        member[i]=j
                        break#jump out of the second rotation, since i has been merged
        #Update member data, to updata some son/grandson points to original father point
        for i in range(Ncd):
            now=int(member[i])
            while now!=member[now]:#If it's not merged points, it will totates to find the father point(merged point)
                now=int(member[now])
            member[i]=now
        return stock,member
