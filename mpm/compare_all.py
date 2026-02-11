#!/usr/bin/env python3
"""
Comprehensive comparison between C++, Python, and Analytical solutions
for MPM 1D vibrating bar problem
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def load_data(filename):
    """Load data from file, handling both formats"""
    data = np.loadtxt(filename, skiprows=1)
    return data

def main():
    print("=" * 80)
    print("MPM 1D: Comprehensive Comparison")
    print("C++ vs Python vs Analytical Solution")
    print("=" * 80)
    
    # Load C++ results
    if not os.path.exists('mpm1D_history.txt'):
        print("ERROR: mpm1D_history.txt not found! Run C++ program first.")
        return
    
    cpp_data = load_data('mpm1D_history.txt')
    time = cpp_data[:, 0]
    x_cpp = cpp_data[:, 1]
    v_cpp = cpp_data[:, 2]
    x_ana = cpp_data[:, 3]
    v_ana = cpp_data[:, 4]
    
    # Load Python results if available
    python_available = os.path.exists('mpm1D_python.txt')
    if python_available:
        py_data = load_data('mpm1D_python.txt')
        x_py = py_data[:, 1]
        v_py = py_data[:, 2]
    else:
        print("\nWarning: mpm1D_python.txt not found. Run Python version first.")
        print("Comparing C++ vs Analytical only.\n")
    
    # Compute errors
    x_err_cpp = np.abs(x_cpp - x_ana)
    v_err_cpp = np.abs(v_cpp - v_ana)
    
    if python_available:
        x_err_py = np.abs(x_py - x_ana)
        v_err_py = np.abs(v_py - v_ana)
        x_diff_cpp_py = np.abs(x_cpp - x_py)
        v_diff_cpp_py = np.abs(v_cpp - v_py)
    
    # Print statistics
    print("\n" + "=" * 80)
    print("ERROR ANALYSIS")
    print("=" * 80)
    
    print("\nC++ Implementation:")
    print(f"  Position - Max: {np.max(x_err_cpp):.6e}, Mean: {np.mean(x_err_cpp):.6e}")
    print(f"  Velocity - Max: {np.max(v_err_cpp):.6e}, Mean: {np.mean(v_err_cpp):.6e}")
    
    if python_available:
        print("\nPython Implementation:")
        print(f"  Position - Max: {np.max(x_err_py):.6e}, Mean: {np.mean(x_err_py):.6e}")
        print(f"  Velocity - Max: {np.max(v_err_py):.6e}, Mean: {np.mean(v_err_py):.6e}")
        
        print("\nC++ vs Python Difference:")
        print(f"  Position - Max: {np.max(x_diff_cpp_py):.6e}, Mean: {np.mean(x_diff_cpp_py):.6e}")
        print(f"  Velocity - Max: {np.max(v_diff_cpp_py):.6e}, Mean: {np.mean(v_diff_cpp_py):.6e}")
    
    print("\n" + "=" * 80)
    print("INITIAL CONDITIONS (t=0)")
    print("=" * 80)
    print(f"Analytical:  x={x_ana[0]:.10f}, v={v_ana[0]:.10f}")
    print(f"C++:         x={x_cpp[0]:.10f}, v={v_cpp[0]:.10f}")
    if python_available:
        print(f"Python:      x={x_py[0]:.10f}, v={v_py[0]:.10f}")
    
    print("\n" + "=" * 80)
    print(f"FINAL STATE (t={time[-1]:.1f})")
    print("=" * 80)
    print(f"Analytical:  x={x_ana[-1]:.10f}, v={v_ana[-1]:.10f}")
    print(f"C++:         x={x_cpp[-1]:.10f}, v={v_cpp[-1]:.10f}, errors: {x_err_cpp[-1]:.6e}, {v_err_cpp[-1]:.6e}")
    if python_available:
        print(f"Python:      x={x_py[-1]:.10f}, v={v_py[-1]:.10f}, errors: {x_err_py[-1]:.6e}, {v_err_py[-1]:.6e}")
    
    # Create comprehensive plots
    n_plots = 4 if python_available else 3
    fig, axes = plt.subplots(n_plots, 2, figsize=(16, 4*n_plots))
    fig.suptitle('MPM 1D: Complete Comparison', fontsize=16, fontweight='bold')
    
    # Plot 1: Position comparison
    ax = axes[0, 0]
    ax.plot(time, x_ana, 'k-', linewidth=2.5, label='Analytical', alpha=0.8)
    ax.plot(time, x_cpp, 'b--', linewidth=1.5, label='C++', alpha=0.8)
    if python_available:
        ax.plot(time, x_py, 'r:', linewidth=1.5, label='Python', alpha=0.8)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Position (m)', fontsize=11)
    ax.set_title('Position vs Time', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Velocity comparison
    ax = axes[0, 1]
    ax.plot(time, v_ana, 'k-', linewidth=2.5, label='Analytical', alpha=0.8)
    ax.plot(time, v_cpp, 'b--', linewidth=1.5, label='C++', alpha=0.8)
    if python_available:
        ax.plot(time, v_py, 'r:', linewidth=1.5, label='Python', alpha=0.8)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Velocity (m/s)', fontsize=11)
    ax.set_title('Velocity vs Time', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Position errors
    ax = axes[1, 0]
    ax.semilogy(time, x_err_cpp, 'b-', linewidth=1.5, label='C++ Error', alpha=0.8)
    if python_available:
        ax.semilogy(time, x_err_py, 'r--', linewidth=1.5, label='Python Error', alpha=0.8)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Absolute Error (m)', fontsize=11)
    ax.set_title('Position Error vs Time', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Velocity errors
    ax = axes[1, 1]
    ax.semilogy(time, v_err_cpp, 'b-', linewidth=1.5, label='C++ Error', alpha=0.8)
    if python_available:
        ax.semilogy(time, v_err_py, 'r--', linewidth=1.5, label='Python Error', alpha=0.8)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Absolute Error (m/s)', fontsize=11)
    ax.set_title('Velocity Error vs Time', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 5: Phase space
    ax = axes[2, 0]
    ax.plot(x_ana, v_ana, 'k-', linewidth=2.5, label='Analytical', alpha=0.8)
    ax.plot(x_cpp, v_cpp, 'b--', linewidth=1.5, label='C++', alpha=0.8)
    if python_available:
        ax.plot(x_py, v_py, 'r:', linewidth=1.5, label='Python', alpha=0.8)
    ax.set_xlabel('Position (m)', fontsize=11)
    ax.set_ylabel('Velocity (m/s)', fontsize=11)
    ax.set_title('Phase Space', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    
    # Plot 6: Relative errors
    ax = axes[2, 1]
    x_rel_err_cpp = x_err_cpp / (np.abs(x_ana) + 1e-10) * 100
    v_rel_err_cpp = v_err_cpp / (np.abs(v_ana) + 1e-10) * 100
    ax.plot(time, x_rel_err_cpp, 'b-', linewidth=1.5, label='C++ Position', alpha=0.8)
    ax.plot(time, v_rel_err_cpp, 'b--', linewidth=1.5, label='C++ Velocity', alpha=0.8)
    if python_available:
        x_rel_err_py = x_err_py / (np.abs(x_ana) + 1e-10) * 100
        v_rel_err_py = v_err_py / (np.abs(v_ana) + 1e-10) * 100
        ax.plot(time, x_rel_err_py, 'r-', linewidth=1.5, label='Python Position', alpha=0.8)
        ax.plot(time, x_rel_err_py, 'r--', linewidth=1.5, label='Python Velocity', alpha=0.8)
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Relative Error (%)', fontsize=11)
    ax.set_title('Relative Errors', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, np.percentile(x_rel_err_cpp, 95) * 1.2])
    
    if python_available:
        # Plot 7: C++ vs Python difference
        ax = axes[3, 0]
        ax.semilogy(time, x_diff_cpp_py, 'g-', linewidth=1.5, label='Position', alpha=0.8)
        ax.semilogy(time, v_diff_cpp_py, 'orange', linewidth=1.5, label='Velocity', alpha=0.8)
        ax.set_xlabel('Time (s)', fontsize=11)
        ax.set_ylabel('Absolute Difference', fontsize=11)
        ax.set_title('C++ vs Python Difference', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Plot 8: Error comparison bar chart
        ax = axes[3, 1]
        categories = ['Pos Max', 'Pos Mean', 'Vel Max', 'Vel Mean']
        cpp_vals = [np.max(x_err_cpp), np.mean(x_err_cpp), 
                    np.max(v_err_cpp), np.mean(v_err_cpp)]
        py_vals = [np.max(x_err_py), np.mean(x_err_py), 
                   np.max(v_err_py), np.mean(v_err_py)]
        
        x_pos = np.arange(len(categories))
        width = 0.35
        ax.bar(x_pos - width/2, cpp_vals, width, label='C++', alpha=0.8)
        ax.bar(x_pos + width/2, py_vals, width, label='Python', alpha=0.8)
        ax.set_ylabel('Error (log scale)', fontsize=11)
        ax.set_title('Error Comparison', fontsize=12, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(categories, rotation=15, ha='right')
        ax.legend(fontsize=10)
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('comparison_complete.png', dpi=300, bbox_inches='tight')
    print(f"\nComprehensive comparison plot saved to: comparison_complete.png")
    
    # Create summary table
    print("\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"{'Metric':<30} {'C++':<15} {'Python':<15} {'Difference':<15}")
    print("-" * 80)
    print(f"{'Position Max Error (m)':<30} {np.max(x_err_cpp):.6e}  ", end="")
    if python_available:
        print(f"{np.max(x_err_py):.6e}  {np.max(x_diff_cpp_py):.6e}")
    else:
        print("N/A             N/A")
    
    print(f"{'Position Mean Error (m)':<30} {np.mean(x_err_cpp):.6e}  ", end="")
    if python_available:
        print(f"{np.mean(x_err_py):.6e}  {np.mean(x_diff_cpp_py):.6e}")
    else:
        print("N/A             N/A")
    
    print(f"{'Velocity Max Error (m/s)':<30} {np.max(v_err_cpp):.6e}  ", end="")
    if python_available:
        print(f"{np.max(v_err_py):.6e}  {np.max(v_diff_cpp_py):.6e}")
    else:
        print("N/A             N/A")
    
    print(f"{'Velocity Mean Error (m/s)':<30} {np.mean(v_err_cpp):.6e}  ", end="")
    if python_available:
        print(f"{np.mean(v_err_py):.6e}  {np.mean(v_diff_cpp_py):.6e}")
    else:
        print("N/A             N/A")
    
    print("=" * 80)
    
    if python_available:
        # Agreement assessment
        max_diff = max(np.max(x_diff_cpp_py), np.max(v_diff_cpp_py))
        if max_diff < 1e-10:
            agreement = "EXCELLENT"
        elif max_diff < 1e-8:
            agreement = "VERY GOOD"
        elif max_diff < 1e-6:
            agreement = "GOOD"
        else:
            agreement = "FAIR"
        
        print(f"\nC++ and Python Implementation Agreement: {agreement}")
        print(f"Maximum difference: {max_diff:.6e}")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
