import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    data = np.loadtxt(file_path)
    return data[:, 0], data[:, 1], data[:, 2]

def calculate_speed_factor(gpu_times, cpu_times):
    return cpu_times / gpu_times

def plot_scaling_test(gpu_times, cpu_times, num_points):
    plt.figure(figsize=(10, 6))
    plt.plot(num_points, gpu_times, linestyle='--', color='blue', label='GPU Time')
    plt.plot(num_points, cpu_times, linestyle='--', color='red',  label='CPU Time')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of Points')
    plt.ylabel('Computational Time (s)')
    plt.title('Scaling Test: GPU vs CPU')
    plt.legend()
    plt.savefig("Scaling_test_GPU.png")

def plot_speed_factor(speed_factor, num_points):
    plt.figure(figsize=(10, 6))
    plt.plot(num_points, speed_factor, linestyle='--', marker='o', color='blue')
    plt.xscale('log')
    plt.xlabel('Number of Points')
    plt.ylabel('Speed Factor (CPU Time / GPU Time)')
    plt.title('Speed Factor: GPU vs CPU')
    plt.savefig('Speed_factor_GPU.png')

# Scaling test plot
file_path = 'GPU_part_Timing.txt'
num_points, gpu_times, cpu_times = read_data(file_path)
plot_scaling_test(gpu_times, cpu_times, num_points)

#Speed factor plot
speed_factor = calculate_speed_factor(gpu_times, cpu_times)
plot_speed_factor(speed_factor, num_points)