import matplotlib.pyplot as plt
import numpy as np

# Read data from file
data = []
with open("interpolation_data_GPU.txt", "r") as file:
    lines = file.readlines()
    for line in lines[1:]:  # Skip the header line
        values = [float(val) for val in line.split()]
        data.append(values)

# Extract data columns
x_data = [row[0] for row in data]
y_data = [row[1] for row in data]
cubic_spline = [row[2] for row in data]

# Read data from file
data_2 = []
with open("interpolation_data_NO_GPU.txt", "r") as file_2:
    lines = file_2.readlines()
    for line in lines[1:]:  # Skip the header line
        values = [float(val) for val in line.split()]
        data_2.append(values)

gsl_spline = [row[3] for row in data_2]

# Create a plot for the original function, cubic spline, and GSL spline
plt.figure(figsize=(12, 6))
# plt.plot(x_data, y_data, label="Original Function", linestyle='--', marker='o', markersize=marker_size)
plt.plot(x_data, cubic_spline, label="Custom Cubic Spline Interpolation", linestyle='-', color='red')
plt.plot(x_data, gsl_spline, label="GSL Spline Interpolation", linestyle='--', color='blue')
plt.xlabel('x_data')
plt.ylabel('Function Values')
plt.title('Comparison of Spline Interpolations')
plt.legend()
plt.savefig("Interpolation_comparison_GPU_10000_data.png")


# Create a plot for the difference between GSL spline and custom interpolation
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
difference = abs(np.array(gsl_spline)) - abs(np.array(cubic_spline))
plt.plot(x_data, difference, label="Difference (GSL - Custom)", linestyle='-', color='green')
plt.xlabel('x_data')
plt.ylabel('GSL Spline - Custom Spline')
plt.legend()


plt.subplot(1, 2, 2)
residuals = (np.array(gsl_spline) - np.array(cubic_spline))/(np.array(gsl_spline))
plt.plot(x_data, residuals, label="Residuals", linestyle='-', color='red')
plt.xlabel('x_data')
plt.ylabel('(GSL Spline - Custom Spline)/(GSL Spline)')
plt.legend()


# Save the figure
plt.tight_layout()  # Adjust layout for better appearance
plt.savefig('Interpolation_comparison_diff_residuals_GPU_10000_data.png')
plt.show()

