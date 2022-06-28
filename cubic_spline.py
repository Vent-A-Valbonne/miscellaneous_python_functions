# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:11:51 2022

@author: Feng_GAO
"""

# import copy
import numpy
import matplotlib.pyplot as plt
# import vav_funs
def cubic_spline_vav(ctrl_p):
    '''
    this work is based on the link：https://bbs.huaweicloud.com/blogs/264151
    
    input is ctrl_p is 2D list-like, in format of ((x_1,y_1),(x_2,y_2),...,(x_n,y_n),(x_n+1,y_n+1)), n+1 points in total
    output is the coefficent array
    the functions are arranaged in the form of n sections of curves:
        
    a1*x**3 + b1*x**2 + c1*x**1 + d1    for 1st section
    a2*x**3 + b2*x**2 + c2*x**1 + d2    for 2nd section
    a3*x**3 + b3*x**2 + c3*x**1 + d3    for 3rd section

    ...
    an*x**3 + bn*x**2 + cn*x**1 + dn    for nth section
    
    '''
    np=len(ctrl_p)
    p_array=numpy.zeros(((np-1)*4,(np-1)*4),dtype=numpy.float32)
    D_array=numpy.zeros(((np-1)*4,),dtype=numpy.float32)
    
    for i in range(0,len(ctrl_p)-1,1):
    # 0 order: the control points, each section has two control points, so totally 2*n functions for the 4*n unknowns
        #first point in section i
        p_array[2*i][4*i+0]=ctrl_p[i][0]**3.0
        p_array[2*i][4*i+1]=ctrl_p[i][0]**2.0
        p_array[2*i][4*i+2]=ctrl_p[i][0]**1.0
        p_array[2*i][4*i+3]=1.0
        #second point in section i
        p_array[2*i+1][4*i+0]=ctrl_p[i+1][0]**3.0
        p_array[2*i+1][4*i+1]=ctrl_p[i+1][0]**2.0
        p_array[2*i+1][4*i+2]=ctrl_p[i+1][0]**1.0
        p_array[2*i+1][4*i+3]=1.0

    for i in range(1,len(ctrl_p)-1,1):
    # 1 order continuity at junction points: totally n+1-2=n-1 functions
        p_array[2*(np-1)+i-1][4*(i-1)+0]=+3.0*ctrl_p[i][0]**2.0
        p_array[2*(np-1)+i-1][4*(i-1)+1]=+2.0*ctrl_p[i][0]**1.0
        p_array[2*(np-1)+i-1][4*(i-1)+2]=+1.0
        p_array[2*(np-1)+i-1][4*(i-1)+4]=-3.0*ctrl_p[i][0]**2.0
        p_array[2*(np-1)+i-1][4*(i-1)+5]=-2.0*ctrl_p[i][0]**1.0
        p_array[2*(np-1)+i-1][4*(i-1)+6]=-1.0
     
    # 2 order continuity at junction points: totally n+1-2=n-1 functions
    for i in range(1,len(ctrl_p)-1,1):
        p_array[3*(np-1)+i-2][4*(i-1)+0]=+6.0*ctrl_p[i][0]**1.0 
        p_array[3*(np-1)+i-2][4*(i-1)+1]=+2.0
        p_array[3*(np-1)+i-2][4*(i-1)+4]=-6.0*ctrl_p[i][0]**1.0 
        p_array[3*(np-1)+i-2][4*(i-1)+5]=-2.0
 
    # boundary condtions at first and last points, 2 functions:
    # here one can change the boundary conditions based on his requirements, refer to the link
    ###### number_1: Natural Spline
    p_array[4*(np-1)-2][0]=6.0*ctrl_p[0][0]
    p_array[4*(np-1)-2][1]=2.0
    p_array[4*(np-1)-1][4*(np-2)]=6.0*ctrl_p[np-1][0]
    p_array[4*(np-1)-1][4*(np-2)+1]=2.0
    # print(p_array)
    
    ###### number_2: Not-a-Knot Spline
    # p_array[4*(np-1)-2][0]=6.0*ctrl_p[1][0]
    # p_array[4*(np-1)-2][4]=-6.0*ctrl_p[1][0]
    # p_array[4*(np-1)-1][4*(np-1)-8]=6.0*ctrl_p[np-2][0]
    # p_array[4*(np-1)-1][4*(np-1)-4]=-6.0*ctrl_p[np-2][0]
    # print(p_array)    
    
    ###### number_3: Quadratic Spline
    # p_array[4*(np-1)-2][0]=1.0e60
    # p_array[4*(np-1)-1][4*(np-1)-4]=1.0e60
    # print(p_array)       
 
    # now the right side
    D_array[0]=ctrl_p[0][1]
    for i in range(0,len(ctrl_p)-2,1):
        D_array[2*i+1]=ctrl_p[i+1][1]
        D_array[2*i+2]=ctrl_p[i+1][1]
    D_array[2*(np-1)-1]=ctrl_p[-1][1]
    # print(D_array)    

    
    answer=numpy.linalg.solve(p_array,D_array)
    return answer


def cubic_curve(x,coe_array,anchor_p):
    '''
    x is the independent scalar for the interpolated curve
    
    coe_array is the (4*n,1) like array for the sectional functions of form:
        a1*x**3 + b1*x**2 + c1*x**1 + d1    for 1st section
        a2*x**3 + b2*x**2 + c2*x**1 + d2    for 2nd section
        a3*x**3 + b3*x**2 + c3*x**1 + d3    for 3rd section
        ...
        an*x**3 + bn*x**2 + cn*x**1 + dn    for nth section
        n sections in total
        
    anchor_p is the anchor points (also call control points) of the curve, in the form of 
    ((x_1,y_1),(x_2,y_2),...,(x_n,y_n),(x_n+1,y_n+1)), n+1 points in total
        
    '''
    x_pos=[anchor_p[i][0] for i in range(len(anchor_p))]
    n=len(x_pos)-1
    for i in range(len(x_pos)-1):
        if x==x_pos[i]:
            y=coe_array[0+i*4]*x**3.0+coe_array[1+i*4]*x**2.0+coe_array[2+i*4]*x+coe_array[3+i*4]
            break
        if (x-x_pos[i])*(x-x_pos[i+1])<0.0:
            y=coe_array[0+i*4]*x**3.0+coe_array[1+i*4]*x**2.0+coe_array[2+i*4]*x+coe_array[3+i*4]
            break
    if x==x_pos[-1]:
        y=coe_array[0+(n-1)*4]*x**3.0+coe_array[1+(n-1)*4]*x**2.0+coe_array[2+(n-1)*4]*x+coe_array[3+(n-1)*4]
 
        
    return y


    


a=[[0.,21.],[1.,24.],[2.,24.],[3.,18],[4.,16.]]

# a=[[-0.8499,0.88],[-0.5245,0.88],[-0.3103,0.88],[-0.16,0.875],[0.,0.86],[0.3103,1.00],[0.4,1.074],[0.5245,1.10],[0.8499,1.10]]
min_x=min([a[i][0] for i in range(len(a))])
max_x=max([a[i][0] for i in range(len(a))])
incre=0.01
x=[incre*i+min_x for i in range(int((max_x-min_x)/incre)+1)]

d=cubic_spline_vav(a)
y=[cubic_curve(i,d,a) for i in x]


plt.figure(figsize = (5,5))#figsize = (5,5) or 'auto'
plt.plot(x,y)
plt.plot([a[i][0] for i in range(len(a))],[a[i][1] for i in range(len(a))],'o')
plt.title('the plotting')
plt.show()





            