import numpy as np
import matplotlib.pyplot as plt


def read_arrays_from_file(filename):
    with open(filename, 'r') as file:
        # Skip the header line
        next(file)
        
        # Read the data and split into columns
        data = [line.split() for line in file]
        data = np.array(data, dtype=float)
        
        # Allocate values into separate arrays
        array1 = data[:, 0]
        array2 = data[:, 1]
        
    return array1, array2

filename = 'out60.txt'  # Replace with the actual input file name or path

try:
    N60, FE60 = read_arrays_from_file(filename)
    print("Arrays read from file successfully:")
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
except ValueError as e:
    print(f"Error: Invalid file format - {str(e)}")

filename = 'out90.txt'  # Replace with the actual input file name or path

try:
    N90, FE90 = read_arrays_from_file(filename)
    print("Arrays read from file successfully:")
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
except ValueError as e:
    print(f"Error: Invalid file format - {str(e)}")

filename = 'out120.txt'  # Replace with the actual input file name or path

try:
    N120, FE120 = read_arrays_from_file(filename)
    print("Arrays read from file successfully:")
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
except ValueError as e:
    print(f"Error: Invalid file format - {str(e)}")

filename = 'out180.txt'  # Replace with the actual input file name or path

try:
    N180, FE180 = read_arrays_from_file(filename)
    print("Arrays read from file successfully:")
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
except ValueError as e:
    print(f"Error: Invalid file format - {str(e)}")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, ax0 = plt.subplots(figsize = (4.0, 4.0), dpi=300)

ax0.plot(N60, FE60, color = 'tab:purple', linestyle = '-', label = r"$60^{\circ}$")#, linestyle = '', marker = 's')
#ax0.plot(N90, FE90, color = 'tab:red', linestyle = '--', label = r"$90^{\circ}$")#, linestyle = '', marker = 's')
#ax0.plot(N120, FE120, color = 'tab:blue', linestyle = ':', label = r"$120^{\circ}$")#, linestyle = '', marker = 's')
#ax0.plot(N180, FE180, color = 'tab:green', linestyle = '-.', label = r"$180^{\circ}$")#, linestyle = '', marker = 's')


ax0.set_xlim([0,100])
ax0.set_ylim([-10,50])
#ax0.set_yticks([])
#ax0.set_xticks([])
ax0.set_xlabel(r"$N_g^V$", fontsize = '14')
ax0.set_ylabel(r"$\Delta \Omega / kT$", fontsize = '14')

plt.legend(fontsize = '14')
plt.tight_layout()
plt.savefig("thetain-out.png")
