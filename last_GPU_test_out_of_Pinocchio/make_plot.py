import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Read data from file
data = []
with open("spline_data.txt", "r") as file:
    lines = file.readlines()
    for line in lines:  # Skip the header line
        values = [float(val) for val in line.split()]
        data.append(values)

# Extract data columns
x_data = [row[0] for row in data]
y_data = [row[1] for row in data]
# cubic_spline = [row[2] for row in data]

# Read data from file
data_2 = []
with open("interpolation_data.txt", "r") as file_2:
    lines = file_2.readlines()
    for line in lines[1:]:  # Skip the header line
        values = [float(val) for val in line.split()]
        data_2.append(values)

actual_values = [row[1] for row in data_2]
evaluated_values =  [row[0] for row in data_2]
cubic_spline = [row[2] for row in data_2]
gsl_spline   = [row[3] for row in data_2]

# Combine the lists into tuples
combined_data = list(zip(evaluated_values, actual_values, cubic_spline, gsl_spline))

# Sort the tuples based on the first element (evaluated_values) in descending order
sorted_data = sorted(combined_data, key=lambda x: x[0], reverse=True)

# Unpack the sorted data into separate lists
sorted_actual_values = [row[1] for row in sorted_data]
sorted_evaluated_values = [row[0] for row in sorted_data]
sorted_cubic_spline = [row[2] for row in sorted_data]
sorted_gsl_spline = [row[3] for row in sorted_data]

# # Create a plot for the original function, cubic spline, and GSL spline
# plt.figure(figsize=(12, 6))
# # plt.plot(x_data, y_data, label="Original Function", linestyle='--', marker='o', markersize=marker_size)
# plt.plot(x_data, y_data, label="Tabulated value", color='black', linestyle='-')
# plt.scatter(sorted_evaluated_values, sorted_gsl_spline, label="GSL Spline Interpolation", marker='^', edgecolors='blue')
# plt.scatter(sorted_evaluated_values, sorted_cubic_spline, label="Custom Spline Interpolation", edgecolors='red', alpha=0.7)
# plt.xlabel('x_data')
# plt.ylabel('Function Values')
# plt.title('Comparison of Spline Interpolations: 30 points')
# plt.legend()
# plt.savefig("Interpolation_comparison_30.png")


# # Create a plot for the difference between GSL spline and custom interpolation
plt.figure(figsize=(12, 6))
# plt.title('Comparison of Spline Interpolations: 30 points')
plt.suptitle('CUSTOM interpolation: 210 points')
plt.subplot(1, 2, 1)
# difference = abs(np.array(sorted_gsl_spline)) - abs(np.array(sorted_cubic_spline))
# plt.scatter(sorted_evaluated_values, difference, label="Difference (GSL - Custom)", linestyle='-', color='green')
# plt.xlabel('x_data')
# plt.ylabel('GSL Spline - Custom Spline')

# # Add average line to the first subplot
# average_difference = np.mean(difference)
# plt.axhline(y=average_difference, color='blue', linestyle='--', label='Average Difference')
# plt.legend()

plt.plot(x_data, y_data, label="Tabulated value", color='black', linestyle='-')
# plt.scatter(sorted_evaluated_values, sorted_actual_values, label="Real values", marker='x', s=80, color='red')
plt.scatter(sorted_evaluated_values, sorted_gsl_spline, label="GSL Spline Interpolation", s=80, edgecolors='blue', facecolors='None')
# plt.scatter(sorted_evaluated_values, sorted_cubic_spline, label="Custom Spline Interpolation", s=80, edgecolors='blue', facecolors='None')
plt.scatter(sorted_evaluated_values, sorted_cubic_spline, label="Custom Spline Interpolation", marker='x', s=80, color='red')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()

plt.subplot(1, 2, 2)
residuals = (np.array(sorted_cubic_spline) - np.array(sorted_gsl_spline))/(np.array(sorted_cubic_spline))
# residuals = (np.array(sorted_gsl_spline) - np.array(sorted_actual_values))/(np.array(sorted_gsl_spline))
# residuals = (np.array(sorted_cubic_spline) - np.array(sorted_actual_values))/(np.array(sorted_cubic_spline))
average_residuals = np.mean(residuals)
plt.scatter(sorted_evaluated_values, residuals, label="Residuals", linestyle='-', color='blue')
formatter = ticker.ScalarFormatter(useMathText=True, useOffset=False)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))  # Adjust if needed
plt.gca().yaxis.set_major_formatter(formatter)
plt.gca().tick_params(axis='y', which='both', labelleft=True)
plt.xlabel('Evaluation point')
plt.ylabel('(CUSTOM Spline -  GSL Spline)/(GSL Spline)')
# plt.ylabel('(GSL Spline -  Actual Value)/(GSL Spline)')
# Add average line to the second subplot
average_residuals = np.mean(residuals)
plt.axhline(y=average_residuals, color='black', linestyle='--', label=f'Average Residuals: {average_residuals*100:.4f} %')
plt.legend()

plt.savefig('210_points_CUSTOM_vs_GSL_yes_boundary.png')

