from ai_recursive import run_ai
import numpy as np
import matplotlib.pyplot as plt
from random import randint


def ComputationComparison():
    naiveAll = []
    naiveMul = []
    aiAll = []
    aiMul = []

    for n in range(2, 100, 1):
        print(n)
        A = generate_random_matrix(n, n, 1, 15)
        B = generate_random_matrix(n, n, 1, 15)

        _, counterAll, counterMul = run_ai(A, B)
        aiAll.append(counterAll)
        aiMul.append(counterMul)

        _, counterAll, counterMul = naive_matrix_multiplication(A, B)
        naiveMul.append(counterMul)
        naiveAll.append(counterAll)
    #
    domain = np.arange(2, 100, 1)
    #plt.plot(domain, naiveAll, color='#8B0000', label='All operands - naive')
    plt.plot(domain, naiveMul, color='#F08080', label='Mul operands - naive')
    #plt.plot(domain, aiAll, color='#006400', label='All operands - AI')
    plt.plot(domain, aiMul, color='#90EE90', label='Mul operands - AI')
    plt.savefig('plot')
    plt.legend()
    plt.show()


def generate_random_matrix(rows, columns, min_val, max_val):
    matrix = np.zeros((rows, columns))
    for a in range(rows):
        for b in range(columns):
            matrix[a][b] = randint(min_val, max_val)

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


if __name__ == '__main__':
    ComputationComparison()
