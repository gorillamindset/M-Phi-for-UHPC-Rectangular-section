import math
import matplotlib.pyplot as plt
### Layer by Layer analysis of Rectangular Beam c/s - UHPC  
### Tensile strength of Concrete is considered as per Bilinear UHPC model
### Tensile stregnth of UHPC after localisation is considered to be Decreasing exponentially
### Units are in N-mm 
Es = 200000             # Young's modulus (Pa)
fc = 150                # Concrete compressive strength (Pa) 
fy = 500                # Yield strength (Pa)
B = 200                 # Width of the beam (m)
D = 500                 # Total depth of the beam (m)
dc = 50                 # Effective cover to the R/f (m)
d = D - dc              # Effective depth of the beam (m)
Ptc = 0/100             # Percentage of compression steel
Ptt = 0.4/100           # Percentage of Tension steel
Asc = Ptc*B*D           # Area of Compression steel reinforcement (m^2)
Ast = Ptt*B*D           # Area of Tension steel reinforcement (m^2)
ec = 0.002                     # Ultimate 44concrete strain
ecmax = 0.0035                 # Ultimate compressive strain in concrete
esmax = 0.12                   # Ultimate strain in steel
fck = (1.25*fc-1.65*5)         # Characteristic comp strength (Pa)
Ec = 0.85*5000*math.sqrt(fc)       # Modulus of elasticity of Concrete (Pa)
fcr = 0.7*math.sqrt(fck)       # modulus of rupture
ecr = fcr/Ec                   # rupture strain (at first tensile crack)
A = B*D                        # Area of cross section
# I = B*D^3/12                 # Moment of Inertia of the section
# Y = D/2                      # Max depth of N.A.
# Z = I/Y                      # Section modulus
Vf = 2                         # Percentage volume of Fibers
Lf = 13                        # Length of Fibers (mm)
Df = 0.2                       # Dia of Fibers (mm)
t = 5                          # Thickness of each layer (mm)
#ectop = 0.002 
n = d/t
print(n)

fcte = 0.65*math.sqrt(fc)
ecte = 1.1*(0.95*fcte)/Ec
RI = Vf*Lf/Df/100
fctl = 0.95*fcte*math.pow((RI),0.25)
ectl = 20*ecte*math.pow((RI),0.25)

def parameters():
    xu = d
    ectop = 0
    list_phi = []
    list_Mu = []
    while ectop < 0.0035:
        while True:
            sum_Cc = 0
            sum_Ct = 0
            sum_Ts = 0
            sum_Mu1 = 0
            counter = 0
            while counter <= n:
                counter += 1
                y_bar = xu - ((counter-1)*t + (t/2))
                ecl = (ectop/xu) * y_bar
                ecl = max(ecl,-ecmax)
                if ecl > 0:
                    fcc = fc * (((2*(ecl/ec))-(math.pow((ecl/ec),2))))
                else:
                    fcc = 0
                if ecl > ecte:
                    fct = fcte/ecte*ecl
                    if ecl > ectl:
                        fct = fcte + ((fctl-fcte)/(ectl-ecte))*(ecl-ecte)
                    else:
                        fct = fctl*(1+(math.pow((ectop-ectl)/ecmax)),3)*math.exp(-2.2*((ectop-ectl)/ecmax))
                                           
                       

                est = ((ectop/xu)*(d-xu))
                if est < (fy/Es):
                    fs = est * Es
                else:
                    fs = fy
                
                Cc = fcc * B * t
                Ct = fct * B * t
                
                Ts = fs * Ast

                Mu1 = (Cc * y_bar) + (Ct * y_bar)
                sum_Mu1 += Mu1
                
                sum_Cc += Cc
                sum_Ct += Ct

            Mu2 = (Ts * (d-xu))
            Mu = (sum_Mu1 + Mu2)/1e6


            if abs(sum_Cc + sum_Ct) - abs(Ts) < 1:
                # print(f"for ectop:{ectop}, parameters: xu:{xu},sum_Cc:{sum_Cc},sum_Ct:{sum_Ct},Mu:{Mu}")

                break
            xu -= 0.001 * d
        ectop += 0.00001
        phi = ectop/xu*(1e6)
        print(f"for ectop:{ectop}, parameters: xu:{xu},sum_Cc:{sum_Cc},sum_Ct:{sum_Ct},Mu:{Mu},phi:{phi}")
        list_phi.append(phi)
        list_Mu.append(Mu)
    plt.plot(list_phi,list_Mu)
    plt.show()
    



parameters()

