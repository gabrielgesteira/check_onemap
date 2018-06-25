# This function was developed to help accounting the recombination frequencies
# outside onemap 'rf_2pts' function, considering ONLY onemap input files
# and the following crosses: RILs ('ril'), back-cross ('bc') and outcross ('outcross'). 
#
# For 'ril' and 'bc' crosses, the function returns a list holding 2 objects:
# the recombination frequency matrix (lower triangular matrix), called 'CC', 
# and the difference between the rf matrix obtained from this function and 
# the 'rf_2pts'(onemap) function, which is called 'DIFF'.
#
# For 'outcross' crosses, the function returns a list holding 8 objetcs:
# 'CC', 'CR', 'RC' and 'RR', which are the recombination frequency matrices 
# (lower triangular matrix) for these four phases. Their respective differences
# when compared to the rf matrices obtained from the 'rf_2pts' (onemap) function
# are called 'DIFFCC', 'DIFFCR', 'DIFFRC' and 'DIFFRR', respectively.

check_onemap = function(file, cross, nround){
  
  # Loading all cross functions (internal)
  cross_ril = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    if (length(t) == 4){
      n1 = t[1,1]
      n2 = t[1,2]
      n3 = t[2,1]
      n4 = t[2,2]
      r = (n2+n3)/(2*(n1+n4))
    } else if (length(t) == 2){
      n1 = t[1]
      n3 = t[2]
      r = n3/(2*n1)
    } else {r=0}
    return(r)
  }
  
  cross_bc = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    if (length(t) == 4){
      n1 = t[1,1]
      n2 = t[1,2]
      n3 = t[2,1]
      n4 = t[2,2]
      r = (n2+n3)/(n1+n2+n3+n4)
    } else if (length(t) == 2){
      n1 = t[1]
      n3 = t[2]
      r = n3/(n1)
    } else {r=0}
    return(r)
  }
  
  ### ML estimators with EM algorithm
  # All crosses and outputs (including linkage phases) follow Maliepaard el at. (1997) and Wu et al. (2002) 
  
  cross_1 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[2,1]
    n4 = t[2,2]
    r1 = (n2 + n3)/n
    r2 = (n1 + n4)/n
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_2 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    r1 = (n3 + n4)/(n1 + n3 + n4 + n6)
    r2 = (n1 + n6)/(n3 + n1 + n6 + n4)
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_3a = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[1,4]
    n5 = t[2,1]
    n6 = t[2,2]
    n7 = t[2,3]
    n8 = t[2,4]
    r1 = (n3 + n4 + n5 + n6)/n
    r2 = (n1 + n2 + n7 + n8)/n
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_3b = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[1,4]
    n5 = t[2,1]
    n6 = t[2,2]
    n7 = t[2,3]
    n8 = t[2,4]
    r1 = (n2 + n4 + n5 + n7)/n
    r2 = (n3 + n1 + n8 + n6)/n
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_4 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[2,1]
    n4 = t[2,2]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((n1*r)/(2-r)) + n2 + ((2*n3*r)/(1+r)))/n
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n3*r)/(2-r)) + n4 + ((2*n1*r)/(1+r)))/n
    }
    
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_5 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    
    r1 = (n2 + n3 + n4)/n
    r2 = (n5 + n6 + n1)/n
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_6 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    
    r1 = (n3 + n5)/(n2 + n3 + n5 + n6)
    r2 = (n2 + n6)/(n3 + n2 + n6 + n5)
    r = matrix(c(r1,r2,r1,r2),2,2)
    return(r)
  }
  
  cross_7 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.001){
      r = r1
      s = 1 - r
      r1 = ((n2 + n4 + n6 + n8 + 2*(n3 + n7) + ((2 * n5 * (r^2)))/(1 - (2 * r * s))))/(2 * n)
    }
    
    # Calculating r2 and r3
    r2 = (0.5 - sqrt((0.25 - ((n1 + n3 + n5 + n7 + n9)/(2 * n)))))
    r3 = r2
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.001){
      r = r4
      s = 1 - r
      r4 = ((n2 + n4 + n6 + n8 + 2*(n1 + n9) + ((2 * n5 * (r^2)))/(1 - (2 * r * s))))/(2 * n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_8 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[1,4]
    n5 = t[2,1]
    n6 = t[2,2]
    n7 = t[2,3]
    n8 = t[2,4]
    n9 = t[3,1]
    n10 = t[3,2]
    n11 = t[3,3]
    n12 = t[3,4]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (n2 + n3 + n5 + n8 + n10 + n11 + 2*(n4 + n9) + (2*(n6 + n7)*(r^2)/(1 - (2*r*s))))/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (n1 + n4 + n6 + n7 + n9 + n12 + 2*(n3 + n10) + (2*(n5 + n8)*(r^2)/(1 - (2*r*s))))/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = (n4 + n1 + n7 + n6 + n12 + n9 + 2*(n2 + n11) + (2*(n8 + n5)*(r^2)/(1 - (2*r*s))))/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (n2 + n3 + n5 + n8 + n10 + n11 + 2*(n1 + n12) + (2*(n6 + n7)*(r^2)/(1 - (2*r*s))))/(2*n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_9 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[2,1]
    n4 = t[2,2]
    n5 = t[3,1]
    n6 = t[3,2]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((2* n1 * r)/(1 + r)) + (2 * n2) + ((n3 * r * (1 + r))/(1 - (r*s))) + n4 + ((2 * n5)/(2 - r)))/(2*n)
    }
    
    # Calculating r2 and r3 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n1 + n5) * r * (1 + r))/(1 - r*s) + n2 + n6 + (2 * n3 * ((1 - s^2))/(1 + 2*r*s)) + ((2 * n4 * r^2)/(1 - (2*r*s))))/(2*n)
    }
    
    r3 = r2
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (((2* n5 * r)/(1 + r)) + (2 * n6) + ((n3 * r * (1 + r))/(1 - (r*s))) + n4 + ((2 * n1)/(2 - r)))/(2*n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_10a = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((n1 + (2 * n4))*r) + n2 + n6 + n8 + (2 * n3) + ((2 * n5 * r^2)/(1 - (2*r*s))) + (n7 *(1 + r)))/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n1 + (2 * n4))*r) + n3 + n5 + n9 + (2 * n2) + ((2 * n6 * r^2)/(1 - (2*r*s))) + (n7 *(1 + r)))/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = (((n7 + (2 * n4))*r) + n9 + n5 + n3 + (2 * n8) + ((2 * n6 * r^2)/(1 - (2*r*s))) + (n1 *(1 + r)))/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (((n7 + (2 * n4))*r) + n2 + n6 + n8 + (2 * n9) + ((2 * n5 * r^2)/(1 - (2*r*s))) + (n1 *(1 + r)))/(2*n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_10b = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((n1 + (2 * n4))*r) + n2 + n6 + n8 + (2 * n3) + ((2 * n5 * r^2)/(1 - (2*r*s))) + (n7 *(1 + r)))/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n7 + (2 * n4))*r) + n9 + n5 + n3 + (2 * n8) + ((2 * n6 * r^2)/(1 - (2*r*s))) + (n1 *(1 + r)))/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = (((n1 + (2 * n4))*r) + n3 + n5 + n9 + (2 * n2) + ((2 * n6 * r^2)/(1 - (2*r*s))) + (n7 *(1 + r)))/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (((n7 + (2 * n4))*r) + n2 + n6 + n8 + (2 * n9) + ((2 * n5 * r^2)/(1 - (2*r*s))) + (n1 *(1 + r)))/(2*n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_11 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[1,4]
    n5 = t[2,1]
    n6 = t[2,2]
    n7 = t[2,3]
    n8 = t[2,4]
    n9 = t[3,1]
    n10 = t[3,2]
    n11 = t[3,3]
    n12 = t[3,4]
    n13 = t[4,1]
    n14 = t[4,2]
    n15 = t[4,3]
    n16 = t[4,4]
    
    r1 = (n2 + n3 + n5 + n8 + n9 + n12 + n14 + n15 + 2*(n4 + n7 + n10 + n13))/(2*n)
    r2 = (n1 + n4 + n6 + n7 + n10 + n11 + n13 + n16 + 2*(n3 + n8 + n9 + n14))/(2*n)
    r3 = (n4 + n1 + n7 + n6 + n11 + n10 + n16 + n13 + 2*(n2 + n5 + n12 + n15))/(2*n)
    r4 = (n2 + n3 + n5 + n8 + n9 + n12 + n14 + n15 + 2*(n1 + n6 + n11 + n16))/(2*n)
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_12 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[2,1]
    n4 = t[2,2]
    n5 = t[3,1]
    n6 = t[3,2]
    n7 = t[4,1]
    n8 = t[4,2]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((2 * n1 * r)/(1 + r)) + (2 * n2) + (((n3 + n5) * r * (1 + r))/(1 - r*s)) + n4 + n6 + ((2 * n7)/(2 - r)))/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((2 * n3 * r)/(1 + r)) + (2 * n4) + (((n1 + n7) * r * (1 + r))/(1 - r*s)) + n2 + n8 + ((2 * n5)/(2 - r)))/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = (((2 * n5 * r)/(1 + r)) + (2 * n6) + (((n7 + n1) * r * (1 + r))/(1 - r*s)) + n8 + n2 + ((2 * n3)/(2 - r)))/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (((2 * n7 * r)/(1 + r)) + (2 * n8) + (((n3 + n5) * r * (1 + r))/(1 - r*s)) + n4 + n6 + ((2 * n1)/(2 - r)))/(2*n)
    }
    
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_13a = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    n10 = t[4,1]
    n11 = t[4,2]
    n12 = t[4,3]
    
    r1 = (n2 + n6 + n7 + n9 + n10 + n11 + 2*(n3 + n5))/(n1 + n4 + n7 + n10 + 2*(n2 + n3 + n5 + n6 + n8 + n9 + n11 + n12))
    r2 = (n3 + n5 + n7 + n8 + n10 + n12 + 2*(n2 + n6))/(n1 + n4 + n7 + n10 + 2*(n3 + n2 + n6 + n5 + n9 + n8 + n12 + n11))
    r3 = (n8 + n12 + n1 + n3 + n4 + n5 + 2*(n9 + n11))/(n7 + n10 + n1 + n4 + 2*(n8 + n9 + n11 + n12 + n2 + n3 + n5 + n6))
    r4 = (n2 + n6 + n4 + n9 + n1 + n11 + 2*(n12 + n8))/(n10 + n7 + n4 + n1 + 2*(n2 + n12 + n8 + n6 + n5 + n9 + n11 + n3))
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_13b = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    n10 = t[4,1]
    n11 = t[4,2]
    n12 = t[4,3]
    
    r1 = (n2 + n6 + n4 + n9 + n10 + n11 + 2*(n3 + n8))/(n1 + n7 + n4 + n10 + 2*(n2 + n3 + n8 + n6 + n5 + n9 + n11 + n12))
    r2 = (n8 + n3 + n7 + n12 + n1 + n5 + 2*(n6 + n11))/(n10 + n4 + n7 + n1 + 2*(n8 + n6 + n11 + n3 + n2 + n12 + n5 + n9))
    r3 = (n3 + n8 + n4 + n5 + n10 + n12 + 2*(n2 + n9))/(n1 + n7 + n4 + n10 + 2*(n3 + n2 + n9 + n8 + n6 + n5 + n12 + n11))
    r4 = (n2 + n6 + n7 + n9 + n1 + n11 + 2*(n12 + n5))/(n10 + n4 + n7 + n1 + 2*(n2 + n12 + n5 + n6 + n8 + n9 + n11 + n3))
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_14 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[2,1]
    n4 = t[2,2]
    
    theta = ((n1 - 2*(n2 + n3) - n4)/(2*n)) + sqrt((((n1 - 2*(n2 + n3) - n4)^2)/((2*n)^2)) + ((2 * n4)/n))
    r1 = 1 - sqrt(theta)
    r2 = 0.5 - sqrt(0.25 - theta)
    r3 = r2
    r4 = sqrt(theta)
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_15a = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((n1 * r * (3 - r))/(2 - r)) + ((n2 * r * (1 + r))/(1 - r*s)) + ((2 * n3)/(2 - r)) + (n4 * (1 + r)) + n5)/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n1 * r * (3 - r))/(2 - r)) + ((n3 * r * (1 + r))/(1 - r*s)) + ((2 * n2)/(2 - r)) + (n4 * (1 + r)) + n5)/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = ((((n1 * r * (3 + r))/(1 + r)) + ((2 * n2 * r)/(1 + r)) + ((n3 * r * (1 + r))/(1 - r*s)) + (n4 * r) + (2 * n5) + n6))/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = ((((n1 * r * (3 + r))/(1 + r)) + ((2 * n3 * r)/(1 + r)) + ((n2 * r * (1 + r))/(1 - r*s)) + (n4 * r) + (2 * n6) + n5))/(2*n)
    }
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_15b = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    
    # Calculating r1 using EM algorithm
    r = 1
    r1 = 0.5
    while (abs(r1 - r) > 0.0001){
      r = r1
      s = 1 - r
      r1 = (((n1 * r * (3 - r))/(2 - r)) + ((n2 * r * (1 + r))/(1 - r*s)) + ((2 * n3)/(2 - r)) + (n4 * (1 + r)) + n5)/(2*n)
    }
    
    # Calculating r2 using EM algorithm
    r = 1
    r2 = 0.5
    while (abs(r2 - r) > 0.0001){
      r = r2
      s = 1 - r
      r2 = (((n1 * r * (3 + r))/(1 + r)) + ((2 * n2 * r)/(1 + r)) + ((n3 * r * (1 + r))/(1 - r*s)) + (n4 * r) + (2 * n5) + n6)/(2*n)
    }
    
    # Calculating r3 using EM algorithm
    r = 1
    r3 = 0.5
    while (abs(r3 - r) > 0.0001){
      r = r3
      s = 1 - r
      r3 = (((n1 * r * (3 - r))/(2 - r)) + ((n3 * r * (1 + r))/(1 - r*s)) + ((2 * n2)/(2 - r)) + (n4 * (1 + r)) + n5)/(2*n)
    }
    
    # Calculating r4 using EM algorithm
    r = 1
    r4 = 0.5
    while (abs(r4 - r) > 0.0001){
      r = r4
      s = 1 - r
      r4 = (((n1 * r * (3 + r))/(1 + r)) + ((2 * n3 * r)/(1 + r)) + ((n2 * r * (1 + r))/(1 - r*s)) + (n4 * r) + (2 * n6) + n5)/(2*n)
    }
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_16a = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    r1 = (n2 + n3 +n4 + n6 + n7 + n8)/(n1 + n2 + n3 + n4 + n7 + 2*(n5 + n6 + n8 + n9))
    r2 = (n2 + n3 +n4 + n5 + n7 + n8)/(n1 + n2 + n3 + n4 + n7 + 2*(n6 + n5 + n8 + n9))
    r3  = (n1 + n5 + n9 + 2*(n6 + n8))/(n1 + n2 + n3 + n4 + n7 + 2*(n5 + n6 + n8 + n9))
    r4  = (n1 + n6 + n8 + 2*(n5 + n9))/(n1 + n2 + n3 + n4 + n7 + 2*(n6 + n5 + n8 + n9))
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_16b = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    r1 = (n2 + n3 +n4 + n6 + n7 + n8)/(n1 + n2 + n3 + n4 + n7 + 2*(n5 + n6 + n8 + n9))
    r2  = (n1 + n5 + n9 + 2*(n6 + n8))/(n1 + n2 + n3 + n4 + n7 + 2*(n5 + n6 + n8 + n9))
    r3 = (n2 + n3 +n4 + n5 + n7 + n8)/(n1 + n2 + n3 + n4 + n7 + 2*(n6 + n5 + n8 + n9))
    r4  = (n1 + n6 + n8 + 2*(n5 + n9))/(n1 + n2 + n3 + n4 + n7 + 2*(n6 + n5 + n9 + n8))
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  cross_17 = function(data){
    data = data[,!((data[1,] == "-") | (data[2,] == "-"))]
    n = ncol(data)
    t = table(as.data.frame(t(data)))
    
    n1 = t[1,1]
    n2 = t[1,2]
    n3 = t[1,3]
    n4 = t[2,1]
    n5 = t[2,2]
    n6 = t[2,3]
    n7 = t[3,1]
    n8 = t[3,2]
    n9 = t[3,3]
    
    r1 = (n3 + n6 + n7 + n8 + (2 * n5))/(n2 + n3 + n4 + n7 + 2*(n5 + n6 + n8 + n9))
    r2 = (n3 + n9 + n4 + n5 + (2 * n8))/(n2 + n3 + n7 + n4 + 2*(n8 + n9 + n5 + n6))
    r3 = (n2 + n5 + n7 + n9 + (2 * n6))/(n3 + n2 + n4 + n7 + 2*(n6 + n5 + n9 + n8))
    r4 = (n2 + n6 + n4 + n8 + (2 * n9))/(n3 + n2 + n7 + n4 + 2*(n9 + n6 + n8 + n5))
    r = matrix(c(r1,r2,r3,r4),2,2)
    return(r)
  }
  
  # Estimating recombination frequency between all markers
  rf_estimate = function(data, cross, nmarkers, nround = 6){
    data_p = as.matrix(data)
    
    if (cross == "outcross"){
      
      # Creating list to hold results
      CC <- matrix(NA,nmarkers,nmarkers)
      rownames(CC) <- colnames(CC) <- substr(as.vector(data_p[,1]),2,4)
      analysis.check <- list()
      analysis.check$'CC' <- CC
      analysis.check$'CR' <- CC
      analysis.check$'RC' <- CC
      analysis.check$'RR' <- CC
      
      # Correcting marker type ID
      data_p[,2] = gsub("\\..*","",data_p[,2])
      
      for(i in 1:(nrow(data_p)-1)){
        for (j in (i+1):nrow(data_p)){
          if (data_p[i,2] == "A" & data_p[j,2] == "A"){
            outp = cross_11(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]    
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "B1"){
            outp = cross_13a(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]    
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "A"){
            outp = cross_13a(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "B2"){
            outp = cross_13b(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "A"){
            outp = cross_13b(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "B3"){
            outp = cross_8(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "A"){
            outp = cross_8(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "C8"){
            outp = cross_12(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "A"){
            outp = cross_12(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "D1"){
            outp = cross_3a(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "A"){
            outp = cross_3a(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "A" & data_p[j,2] == "D2"){
            outp = cross_3b(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "A"){
            outp = cross_3b(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "B1"){
            outp = cross_16a(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "B2"){
            outp = cross_17(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "B1"){
            outp = cross_17(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "B3"){
            outp = cross_10a(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "B1"){
            outp = cross_10a(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "C8"){
            outp = cross_15a(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "B1"){
            outp = cross_15a(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "D1"){
            outp = cross_5(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "B1"){
            outp = cross_5(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B1" & data_p[j,2] == "D2"){
            outp = cross_6(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "B1"){
            outp = cross_6(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "B2"){
            outp = cross_16b(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "B3"){
            outp = cross_10b(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "B2"){
            outp = cross_10b(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "C8"){
            outp = cross_15b(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "B2"){
            outp = cross_15b(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "D1"){
            outp = cross_6(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "B2"){
            outp = cross_6(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B2" & data_p[j,2] == "D2"){
            outp = cross_5(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "B2"){
            outp = cross_5(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "B3"){
            outp = cross_7(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "C8"){
            outp = cross_9(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "B3"){
            outp = cross_9(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "D1"){
            outp = cross_2(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "B3"){
            outp = cross_2(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "B3" & data_p[j,2] == "D2"){
            outp = cross_2(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "B3"){
            outp = cross_2(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "C8"){
            outp = cross_14(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "D1"){
            outp = cross_4(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "C8"){
            outp = cross_4(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "C8" & data_p[j,2] == "D2"){
            outp = cross_4(data_p[c(j,i),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "C8"){
            outp = cross_4(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "D1"){
            outp = cross_1(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
          else if (data_p[i,2] == "D1" & data_p[j,2] == "D2"){
            outp = matrix(c(0,0,0,0),2,2)
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
            warning("Missing cross information", call. = F)
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "D1"){
            outp = matrix(c(0,0,0,0),2,2)
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
            warning("Missing cross information", call. = F)
          }
          else if (data_p[i,2] == "D2" & data_p[j,2] == "D2"){
            outp = cross_1(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp[1,1]
            analysis.check[[2]][i,j] <- analysis.check[[2]][j,i] <- outp[2,1]
            analysis.check[[3]][i,j] <- analysis.check[[3]][j,i] <- outp[1,2]
            analysis.check[[4]][i,j] <- analysis.check[[4]][j,i] <- outp[2,2]
          }
        }
      }
      analysis.check[] <- lapply(analysis.check,round,nround)
      return(analysis.check)
      
    } else {
      # Creating list to hold results
      CC <- matrix(NA,nmarkers,nmarkers) 
      rownames(CC) <- colnames(CC) <- data[,1] # substr(as.vector(data_p[,1]),2,4)
      analysis.check <- list()
      analysis.check$'CC' <- CC
      
      if (cross == "ril"){
        for(i in 1:(nrow(data_p)-1)){
          for (j in (i+1):nrow(data_p)){
            outp = cross_ril(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp
          }
        }
      } else if (cross == "bc"){
        for(i in 1:(nrow(data_p)-1)){
          for (j in (i+1):nrow(data_p)){
            outp = cross_bc(data_p[c(i,j),-c(1,2)])
            analysis.check[[1]][i,j] <- analysis.check[[1]][j,i] <- outp
          }
        }
      }
      analysis.check[] <- lapply(analysis.check,round,nround)
      return(analysis.check)
    }
  }
  
  # Reading .raw file
  data = read.table(file, skip = 3) # Only considering Onemap files
  data = as.matrix(data)
  nmarkers = nrow(data)
  
  # Estimating recombination frequency outside onemap
  check_twopts = rf_estimate(data, cross, nmarkers, nround)
  
  if (cross == "outcross"){
    # Correcting missing values
    check_twopts$CC[which(is.na(check_twopts$CC))] <- 0
    check_twopts$CR[which(is.na(check_twopts$CR))] <- 0
    check_twopts$RC[which(is.na(check_twopts$RC))] <- 0
    check_twopts$RR[which(is.na(check_twopts$RR))] <- 0
    
    # Proceeding analysis in onemap
    onemap_data = onemap::read_onemap(inputfile = file)
    twopts = onemap::rf_2pts(onemap_data, LOD = onemap::suggest_lod(onemap_data), max.rf = 0.5)
    
    # Verifying differences
    check_twopts$'DIFFCC' <- round(twopts$analysis$CC - check_twopts$CC, nround)
    check_twopts$'DIFFCR' <- round(twopts$analysis$CR - check_twopts$CR, nround)
    check_twopts$'DIFFRC' <- round(twopts$analysis$RC - check_twopts$RC, nround)
    check_twopts$'DIFFRR' <- round(twopts$analysis$RR - check_twopts$RR, nround)
    
    cat("Done! Obs: D1 x D2 crosses generate warning messages.")
  } else {
    # Correcting rf missing values and those bigger than 0.5
    check_twopts$CC[which(check_twopts$CC > 0.5)] <- 0.5
    check_twopts$CC[which(is.na(check_twopts$CC))] <- 0
    
    # Proceeding analysis in onemap
    onemap_data = onemap::read_onemap(inputfile = file)
    twopts = onemap::rf_2pts(onemap_data, LOD = onemap::suggest_lod(onemap_data), max.rf = 0.5)
    
    # Verifying differences
    check_twopts$'DIFF' <- round(twopts$analysis - check_twopts$CC, nround)
    cat("Done!")
  }
  return(check_twopts)
}