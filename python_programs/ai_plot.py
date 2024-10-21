import time
import numpy as np
import matplotlib.pyplot as plt
from random import randint
from ai_recursive import run_ai


def plot_time():
    naive_times = []
    ai_times = []

    for n in range(2, 200, 1):
        print(f'Processing matrix size {n}x{n}...')
        A = generate_random_matrix(n, n, 0.00000001, 1.0)
        B = generate_random_matrix(n, n, 0.00000001, 1.0)

        # Measure time for AI-based multiplication
        start_time = time.time()
        run_ai(A, B)  # Run the AI multiplication
        ai_time = time.time() - start_time
        ai_times.append(ai_time)

        # Measure time for naive multiplication
        start_time = time.time()
        naive_matrix_multiplication(A, B)  # Run the naive multiplication
        naive_time = time.time() - start_time
        naive_times.append(naive_time)

    domain = np.arange(2, 200, 1)

    # Plot the times
    plt.plot(domain, naive_times, color='purple', linestyle='--', marker='o', label='Naive Multiplication Time')
    plt.plot(domain, ai_times, color='forestgreen', linestyle='-', marker='^', label='AI Multiplication Time')

    plt.legend(title='Execution Time', loc='upper left', fancybox=True, shadow=True, fontsize=10)
    plt.title('Comparison of Execution Times')
    plt.xlabel('Matrix Size (n)')
    plt.ylabel('Time (seconds)')
    plt.grid(True)

    plt.savefig('execution_time_comparison.png')
    plt.show()


def plot_flops():
    naive_all = []
    naive_mul = []
    ai_all = []
    ai_mul = []

    for n in range(2, 50, 1):
        print(n)
        A = generate_random_matrix(n, n, 0.00000001, 1.0)
        B = generate_random_matrix(n, n, 0.00000001, 1.0)

        _, counter_all, counter_mul = run_ai(A, B)
        ai_all.append(counter_all)
        ai_mul.append(counter_mul)

        _, counter_all, counter_mul = naive_matrix_multiplication(A, B)
        naive_mul.append(counter_mul)
        naive_all.append(counter_all)

    domain = np.arange(2, 50, 1)
    plt.plot(domain, naive_all, color='purple', linestyle='--', marker='o',
             label='All operands (Naive)')
    plt.plot(domain, naive_mul, color='skyblue', linestyle='-', marker='s',
             label='Multiplication operands (Naive)')
    plt.plot(domain, ai_all, color='forestgreen', linestyle='-.', marker='^',
             label='All operands (AI)')
    plt.plot(domain, ai_mul, color='coral', linestyle=':', marker='x',
             label='Multiplication operands (AI)')

    plt.legend(title='Computation Types', loc='upper left', fancybox=True,
               shadow=True,fontsize=10)

    plt.title('Comparison of Computation Operations')
    plt.xlabel('Matrix Size (n)')
    plt.ylabel('Number of Operations')
    plt.grid(True)

    plt.savefig('computation_comparison.png')
    plt.show()


def generate_random_matrix(rows, columns, min_val, max_val):
    matrix = np.random.uniform(low=min_val, high=max_val, size=(rows, columns))

    return np.array(matrix)


def naive_matrix_multiplication(A, B):
    counter_all = 0
    counter_mul = 0

    if len(A[0]) != len(B):
        raise ValueError("Matrices size mismatch")

    result = np.zeros((len(A), len(B[0])))

    for a in range(len(A)):
        for b in range(len(B[0])):
            for c in range(len(B)):
                counter_all += 2
                counter_mul += 1
                result[a][b] += A[a][c] * B[c][b]

    return result, counter_all, counter_mul


def estimate():
    avg = 0
    cnt = 0
    for n in range(20, 50, 1):
        cnt += 1
        print(n)
        A = generate_random_matrix(n, n, 0.00000001, 1.0)
        B = generate_random_matrix(n, n, 0.00000001, 1.0)

        start_time = time.time()
        run_ai(A, B)
        ai_time_n = time.time() - start_time

        A = generate_random_matrix(2*n,2*n, 0.00000001, 1.0)
        B = generate_random_matrix(2*n, 2*n, 0.00000001, 1.0)
        start_time = time.time()
        run_ai(A, B)
        ai_time_2n = time.time() - start_time
        avg += ai_time_2n / ai_time_n

    print(avg / cnt)


if __name__ == '__main__':
    estimate()