def get_quad(ang, x, y):
    from numpy import pi
    #ang = bs(-pi*0.7)%(2*pi) * sign(-pi*0.7)
    if (ang > 0 and ang <= pi/2):  #>0 to 90: Q1
        return(x, y)
    elif (ang > pi/2 and ang<= pi): #>90 to 180: Q2
        return(-x, y)
    elif (ang > pi and ang<= 3*pi/2): #>180 to 270: Q3
        return(-x, -y)
    elif (ang > 3*pi/2 and ang<= 2*pi): #>270 to 360: Q4
        return(x,-y)
    elif (ang <= 0 and ang > -pi/2):   #>-90 to 0: Q4
        return(x,-y)
    elif (ang <= -pi/2 and ang > -pi):   #>-180 to -90: Q3
        return(-x,-y)
    elif (ang <= -pi and ang > -3*pi/2):   #>-270 to -180 : Q2
        return(-x,y)
    elif (ang <= -3*pi/2 and ang > -2*pi):   #>-360 to -270: Q1
        return(x,y)
