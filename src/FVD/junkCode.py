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


def fixCoeffsP (coeff,side):
    i = np.size(coeff, 0)
    j = np.size(coeff, 1)

    if side in 'W':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  != (i-1) and n == 0 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == 0 :
                    coeff[m][n] = 0.0

    if side in 'E':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  != (i-1) and n == j-1 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == j-1 :
                    coeff[m][n] = 0.0

    if side in 'S':
        for m in range(i):  # loop through rows
            for n in range(j):  # loop through columns
                if m  == (i-1) and n != j-1 :
                    coeff[m][n] = 0.0

                if m  == (i-1) and n == j-1 :
                    coeff[m][n] = 0.0
    return coeff


def cloneRowMat(F, G, l, p):
    """A function to clone pth row from a matrix G to into the lth row of a matrix F (p > l)"""
    j = np.size(F, 1)
    for m in range(p):  # loop through rows
        for n in range(j):  # loop through columns
            if (m == p - 1):
                F[m][n] = G[l - 1][n]
    return F


def cloneColMat(F, G, l, p):
    """A function to clone lth column into the pth one (p > l)"""
    Fnew = F
    j = np.size(Fnew, 0)
    for m in range(j):  # loop through rows
        for n in range(p):  # loop through columns
            if (n == p - 1):
                Fnew[m][n] = G[m][l - 1]
    return Fnew


if (m > 0 and m < i - 1 and n == 1):  # Boundary face (WEST):  # first grid nodes
    PP = Px[m][n]
    PW = Px[m][n - 1]
    PE = Px[m][n + 1]
    PEE = Px[m][n + 1]

    # East faces
    coeff1e = rhicE(PEE, PE, PP, PW)
    coeff2e = (dy[m][n] ** 2) * rho / (4.0 * aEe[m][n])
    pcorre[m][n] = coeff1e / coeff2e

if (m == i - 2 and n > 0 and n < j - 1):  # Boundary face (SOUTH):  # first grid nodes
    PP = Px[m][n]
    PN = Px[m - 1][n]
    PNN = Px[m - 1][n]
    PS = Px[m + 1][n]

    # North faces
    coeff1n = rhicN(PNN, PN, PP, PS)  # Pn = Pp (zero gradient bc)
    coeff2n = (dx[m][n] ** 2) * rho / (4.0 * aNn[m][n])
    pcorrn[m][n] = coeff1n / coeff2n


for i in range(1, nx-1):
       for j in range(1, ny-1):
           tmp = u[i,j]
           u[i,j] = ((u[i-1, j] + u[i+1, j])*dy2 +
                    (u[i, j-1] + u[i, j+1])*dx2)*dnr_inv
           diff = u[i,j] - tmp
                   err += diff*diff


# The actual iteration
           u[1:-1, 1:-1] = ((u[0:-2, 1:-1] + u[2:, 1:-1])*dy2 +
                            (u[1:-1,0:-2] + u[1:-1, 2:])*dx2)*dnr_inv




#
for m in range(i):  # loop through rows
    for n in range(j):  # loop through columns
        if (m != 0 and n != 0 and m != (i - 1) and n != (j - 1)):  # Internal nodes:

            uN = u[m - 1][n]
            uS = u[m + 1][n]
            uE = u[m][n + 1]
            uW = u[m][n - 1]

            u[1:-1, 1:-1] = (u[1:-1,0:-2] * aW[1:-1, 1:-1] + u[1:-1, 2:] * aE[1:-1, 1:-1] +  u[2:, 1:-1] * aS[1:-1, 1:-1] + \
                             u[0:-2, 1:-1] * aN[1:-1, 1:-1] + SUx[1:-1, 1:-1]) / aP[1:-1, 1:-1]
            sumRes = sumRes + (uW * aW[m][n] + uE * aE[m][n] + uS * aS[m][n] + \
                               uN * aN[m][n] + SUx[m][n]) - u[m][n] * aP[m][n]