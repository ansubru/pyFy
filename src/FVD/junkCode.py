def DDChk(i,j,A):
    sum = np.zeros((i, 1), dtype=float)
    flag = []
    for m in range(i):
        for n in range(j):
            if m != n:
                sum[m] = sum[m] + abs(A[m][n])

    for m in range(i):
        for n in range(j):
            if m == n:
                if A[m][n] > sum[m]:
                    flag.append("DD")
                else:
                    flag.append("NXD")

    for m in range(flag.__len__()):
        if flag[m] in "NXD":
            print ("Critical error! Matrix is not diagonally dominant at row ",m)
            state = "fail"
        else:
            state = "pass"
    print state
    return state

def gaussSeidelP(self, matP, aWinp, aEinp, aNinp, aSinp, aPinp, Bin , iterinp):
        """Solves gauss seidel for a fixed number of iterations."""
        print("Trying to solve equations for P using the Gauss-Seidel2 method")

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")

        #Initialize all relevant variables:u, v, p etc.
        u = 1.0*matP
        usolve, uP, uW, uE, uN, uS = 0.0*matP, 0.0*matP, 0.0*matP, 0.0*matP,0.0*matP, 0.0*matP
        aW = np.array(aWinp)
        aE = np.array(aEinp)
        aS = np.array(aSinp)
        aN = np.array(aNinp)
        aP = np.array(aPinp)
        b = np.array(Bin)
        i = np.size(matP, 0)
        j = np.size(matP, 1)

        iter = 0 # Mock up error parameter only to run the While

        while iter > iterinp:
            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == 0 and n == 0):  # cells in the top left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]


                    elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]

                    elif m == 0 and n == (j-1):  # Top right edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = u[m][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == (i-1) and n == 0:  # bottom left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m][n]
                        uN[m][n] = u[m-1][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n]
                       uS[m][n] = u[m][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]


                    elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + b[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

            iter =+ 1
        return usolve

def gaussSeidel2v(self, matV, aWinp, aEinp, aNinp, aSinp, aPinp, SUyin , iterinp):
        """Solves gauss seidel for a fixed number of iterations."""
        print("Trying to solve V velocities ...")

        ###############------------CREATE IO OBJECT--------------################
        from IO import IO
        IO_obj = IO("random")
        #Boundary conditions
        VA = IO_obj.VA
        VB = IO_obj.VB
        VC = IO_obj.VC
        VD = IO_obj.VD

        #Initialize all relevant variables:u, v, p etc.
        u = 1.0*matV
        usolve, uP, uW, uE, uN, uS = 0.0*matV, 0.0*matV, 0.0*matV, 0.0*matV,0.0*matV, 0.0*matV
        aW = np.array(aWinp)
        aE = np.array(aEinp)
        aS = np.array(aSinp)
        aN = np.array(aNinp)
        aP = np.array(aPinp)
        SUy = np.array(SUyin)
        i = np.size(matV, 0)
        j = np.size(matV, 1)

        iter = 0 # Mock up error parameter only to run the While

        while iter > iterinp:
            for m in range(i):  # loop through rows
                for n in range(j):  # loop through columns
                    if (m == 0 and n == 0):  # cells in the top left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = VA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]

                    elif m == 0 and n != (j-1):  # First row bordering the BC --> B (sans edges)

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]

                    elif m == 0 and n == (j-1):  # Top right edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = u[m][n-1]
                        uE[m][n] = VC
                        uS[m][n] = u[m+1][n]
                        uN[m][n] = VB

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif m == (i-1) and n == 0:  # bottom left edge

                        uP[m][n] = u[m][n]
                        uW[m][n] = VA
                        uE[m][n] = u[m][n+1]
                        uS[m][n] = VD
                        uN[m][n] = u[m-1][n]

                        #Solve gauss seidel 2
                        usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                        u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n != (j-1)):  # Bottom row bordering the BC --> D  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = VD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m != (i-1) and n == 0):  # First column bordering the BC --> A  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = VA
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif (m == (i-1) and n == (j-1)):  # Bottom right edge

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = VC
                       uS[m][n] = VD
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]


                    elif (m != (i-1) and n == (j-1)):  # Right (last) column bordering the BC --> C  (sans edges)

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = VC
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

                    elif m != 0 and n != 0 and m != (i-1) and n != (j-1):  # All other elements

                       uP[m][n] = u[m][n]
                       uW[m][n] = u[m][n-1]
                       uE[m][n] = u[m][n+1]
                       uS[m][n] = u[m+1][n]
                       uN[m][n] = u[m-1][n]

                       #Solve gauss seidel 2
                       usolve[m][n] = (uW[m][n]*aW[m][n] + uE[m][n]*aE[m][n] + uS[m][n]*aS[m][n] + uN[m][n]*aN[m][n] + SUy[m][n])/aP[m][n]
                       u[m][n] = usolve[m][n]

            iter =+ 1

        return usolve

#Create array with x,y co-ordinates from generated grid
        i = np.size(grid_x,0)
        j = np.size(grid_y,1)

        for m in range(0,j):
            for n in range(0,i):
                if(m==n==0):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = y
                elif(m == 0 and n < i):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = y
                elif(m > 0 and n > 0):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = (y - ((m)*delta_y))
                elif(n == 0 and m == j):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = 0.0
                elif(n == 0 and m < j):
                    grid_x[m][n] = 0.0
                    grid_y[m][n] = (y - ((m)*delta_y))
                elif(n == i and m == j):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = 0.0
                elif(n == i and m == 0):
                    grid_x[m][n] = (n)*delta_x
                    grid_y[m][n] = y

