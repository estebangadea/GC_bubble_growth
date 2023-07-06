#Free energy of a bubble.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, log, acos, sin, asin, exp

def read_file(filename):
    variables = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                name, value = line.split('=')
                name = name.strip()
                value = value.strip()
                variables[name] = float(value)
    return variables

def write_arrays_to_file(filename, array1, array2):
    if array1.shape != array2.shape:
        raise ValueError("Arrays must have the same shape.")
    if len(array1.shape) != 1 or len(array2.shape) != 1:
        raise ValueError("Arrays must be 1-dimensional.")

    with open(filename, 'w') as file:
        file.write("N  FE\n")
        for n, fe in zip(array1, array2):
            file.write(f"{n}  {fe}\n")

def get_N(V, T, P): #Rimbach, H., Chatterjee, N.D. Equations of state for H2, H2O, and H2-H2O fluid mixtures at temperatures above 0.01° C and at high pressures. Phys Chem Minerals 14, 560–569 (1987). https://doi.org/10.1007/BF00308292

	a = 0.2476*100000/1000**2
	b = 0.02661/1000

	R = 8.31446261815324

	C03 = - a * b / V**2
	C02 = a / V
	C01 = -1 * (R * T + P * b)
	C00 = P * V
	roots = np.roots([C03, C02, C01, C00])
	#print(roots)
	Vscaled = V / np.real(roots[2]) / b
	Tscaled = T * R * b / a 
	Pscaled = P * b * b / a
	lnfug = 1 / (Vscaled - 1) - 2 / (Tscaled * Vscaled) + log(Tscaled / (Pscaled * (Vscaled -1)))
	fug = exp(lnfug) * P
	return 6.02e23*np.real(roots[2]), fug

def get_geom(P, gamma0, Pext, theta_in, amax, R):
	gamma = gamma0
	r = (R + 2*gamma / (P-Pext))/2
	a = sin(theta_in) * r
	theta = theta_in
	h = np.roots([1, -2*r, a**2])[1]
	if a>amax:
		a = amax
		h = np.roots([1, -2*r, a**2])[1]
		if h<=r:
			theta = asin(amax/r)
		elif h>r:
			theta = np.pi - asin(amax/r)		
	h = np.roots([1, -2*r, a**2])[1]
	V = 1.0/6.0 * np.pi * h * (3*a**2 + h**2)
	area = 2*np.pi*R**2*(1-cos(theta))
	base = np.pi * a**2

	return r, theta, a, h, V, area, base



def free_energy(N, R, theta, vol, area, base, gamma0, kh, pext, concentration, fug):
	gamma = gamma0
	c1 = N * 8.314462618*298/6.02e23 * log(kh * (fug)/concentration)
	c2 = area * gamma
	c3 = -1 * base * cos(theta_in) * gamma
	c4 = - 2 * gamma / R * vol
	return c1, c2, c3, c4

def free_energy2(N, R, theta, vol, area, base_in, base, gamma0, kh, pext, concentration, fug):
	gamma = gamma0
	c1 = N * 8.314462618*298/6.02e23 * log(kh * (fug)/concentration)
	c2 = area * gamma
	c3 = -1 * base_in * cos(theta_in) * gamma - 1 * (base-base_in) * cos(theta_out) * gamma
	#c3 = -1 * base * cos(theta_out) * gamma
	c4 = - 2 * gamma / R * vol
	return c1, c2, c3, c4

def get_geom_ideal(N, theta_in, pext, gamma, maxr):
	RT = 8.314462618*298
	a = pext * np.pi/3 * (cos(theta_in)-1)**2*(2+cos(theta_in))/RT
	b = 2 * gamma * np.pi/3 * (cos(theta_in)-1)**2*(2+cos(theta_in))/RT
	c = - N/6.02e23
	R = np.real(np.roots([a, b, 0, c])[2])
	baser = sin(theta_in) * R
	if baser <= maxr:
		theta = theta_in
		h = np.roots([1, -2*R, baser**2])[1]
	else:
		c5 = pext
		c4 = 4 * gamma
		c3 = 4 * pext * maxr**2
		c2 = 12 * gamma * maxr**2 - 6 * N/6.02e23 * RT / np.pi
		c1 = 3 * pext * maxr**4
		c0 = -6 * N/6.02e23 * RT * maxr**2 / np.pi
		for root in np.roots([c5,c4,c3,c2,c1,c0]):
			if np.real(root) > 0:
				h = np.real(root)
		#print(np.roots([c5,c4,c3,c2,c1,c0]))
		R = (maxr**2 + h**2) / (2*h)
		if h<=R:
			theta = asin(maxr/R)
		elif h>R:
			theta = np.pi - asin(maxr/R)
		baser = maxr

	vol = np.pi/3 * R**3 * (cos(theta)-1)**2*(2+cos(theta))
	area = 2*np.pi*R**2*(1-cos(theta))
	base = np.pi * baser**2
	#print(theta)
	return R, area, base, vol, theta, baser, h

def try_c(concentration):

	a = []
	h = []
	r = []
	theta = []
	V = []
	A = []
	B = []
	P = []
	N = []
	FE = []
	FEc1 = []
	FEc2 = []
	FEc3 = []
	FEc4 = []


	for ii in range(160):
		ai = (ii+40) * maxr/200.0
		a.append(ai)
		ri = ai/sin(theta_in)
		r.append(ri)
		hi = ri - ri * cos(theta_in)
		h.append(hi)
		thetai = theta_in
		theta.append(thetai)
		Vi = np.pi/3.0 * ri**3 * (2 + cos(theta_in)) * (1-cos(theta_in))**2
		V.append(Vi)
		Ai = 2 * np.pi * ri**2 * (1-cos(theta_in))
		A.append(Ai)
		Bi = np.pi * ai**2
		B.append(Bi)
		Pi = (pext + 2*gamma/ri) / (1 + 2.4673e-10 / ri)
		P.append(Pi)
		Ni, fugi = get_N(Vi, 298.0, Pi)
		N.append(Ni)
		FEi = free_energy(Ni, ri, theta_in, Vi, Ai, Bi, gamma, kh, pext, concentration, fugi)
		FEc1.append(FEi[0]/(1.380649e-23*298))
		FEc2.append(FEi[1]/(1.380649e-23*298))
		FEc3.append(FEi[2]/(1.380649e-23*298))
		FEc4.append(FEi[3]/(1.380649e-23*298))
		FE.append(np.sum(FEi)/(1.380649e-23*298))
	
	ii=0
	while thetai < theta_out:
		ii+=1
		hi = h[159] + (ii+1) * h[159]/200.0
		h.append(hi)
		ai = maxr
		a.append(ai)
		ri = (ai**2 + hi**2)/(2*hi)
		r.append(ri)
		if hi<=ri:
			thetai = asin(ai/ri)
		elif hi>ri:
			thetai = np.pi - asin(ai/ri)
		theta.append(thetai)
		Vi = np.pi/3.0 * ri**3 * (2 + cos(thetai)) * (1-cos(thetai))**2
		V.append(Vi)
		Ai = 2 * np.pi * ri**2 * (1-cos(thetai))
		A.append(Ai)
		Bi = np.pi * ai**2
		B.append(Bi)
		Pi = (pext + 2*gamma/ri) / (1 + 2.4673e-10 / ri)
		P.append(Pi)
		Ni, fugi = get_N(Vi, 298.0, Pi)
		N.append(Ni)
		FEi = free_energy(Ni, ri, thetai, Vi, Ai, Bi, gamma, kh, pext, concentration, fugi)
		FEc1.append(FEi[0]/(1.380649e-23*298))
		FEc2.append(FEi[1]/(1.380649e-23*298))
		FEc3.append(FEi[2]/(1.380649e-23*298))
		FEc4.append(FEi[3]/(1.380649e-23*298))
		FE.append(np.sum(FEi)/(1.380649e-23*298))


	for ii in range(300):
		ai = maxr + ii * maxr/200.0
		a.append(ai)
		ri = ai/sin(theta_out)
		#print(ri)
		r.append(ri)
		hi = ri - ri * cos(theta_out)
		h.append(hi)
		thetai = theta_out
		theta.append(thetai)
		Vi = np.pi/3.0 * ri**3 * (2 + cos(thetai)) * (1-cos(thetai))**2
		V.append(Vi)
		Ai = 2 * np.pi * ri**2 * (1-cos(thetai))
		A.append(Ai)
		Bi_in = np.pi * maxr**2
		Bi = np.pi * ai**2
		B.append(Bi)
		Pi = (pext + 2*gamma/ri) / (1 + 2.4673e-10 / ri)
		P.append(Pi)
		Ni, fugi = get_N(Vi, 298.0, Pi)
		N.append(Ni)
		FEi = free_energy2(Ni, ri, thetai, Vi, Ai, Bi_in, Bi, gamma, kh, pext, concentration, fugi)
		FEc1.append(FEi[0]/(1.380649e-23*298))
		FEc2.append(FEi[1]/(1.380649e-23*298))
		FEc3.append(FEi[2]/(1.380649e-23*298))
		FEc4.append(FEi[3]/(1.380649e-23*298))
		FE.append(np.sum(FEi)/(1.380649e-23*298))

	return np.array(N), np.array(FE), np.array(FEc1), np.array(FEc2), np.array(FEc3), np.array(FEc4)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, ax0 = plt.subplots(figsize = (4.0, 4.0), dpi=300)


filename = 'input.dat'  # Replace with your file name or path

try:
    values = read_file(filename)
    maxr = values['maxr']
    theta_in = np.pi*values['theta_in']/180.0
    theta_out = np.pi*values['theta_out']/180.0
    concentration = values['concentration']
    kh = values['kh']
    pext = values['pext']
    gamma = values['gamma']

    # Now you can use the variables as needed
    print(maxr)
    print(theta_in)
    print(theta_out)
    print(concentration)
    print(kh)
    print(pext)
    print(gamma)
except FileNotFoundError:
    print(f"File '{filename}' not found.")
except ValueError:
    print("Error: Invalid file format.")

N, FE, FEc1, FEc2, FEc3, FEc4 = try_c(concentration)

for jj in range(len(FE)-2):
	if FE[jj]< FE[jj+1] and FE[jj+1]>FE[jj+2]:
		print(("N: %i  FE: %f")%(N[jj+1],FE[jj+1]))
		break

filename = 'out.txt'  # Replace with the desired output file name or path

try:
    write_arrays_to_file(filename, N, FE)
    print(f"Arrays written to '{filename}' successfully.")
except ValueError as e:
    print(f"Error: {str(e)}")


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2, figsize = (6.0, 6.0), dpi=300)


ax0.plot(N, FEc1)#, linestyle = '', marker = 's')
ax1.plot(N, FEc2)#, linestyle = '', marker = 's')
ax2.plot(N, FEc3)#, linestyle = '', marker = 's')
ax3.plot(N, FEc4)#, linestyle = '', marker = 's')
#ax0.set_xticks([])

#plt.legend(fontsize = '18')
plt.tight_layout()
plt.savefig("asdasd.png")