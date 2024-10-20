import unittest
import numpy as np
from random import randint
from ai_recursive import run_ai


class AITest(unittest.TestCase):
    def test_MultiplicationTest(self):
        for n in range(2, 30):
            A = generate_random_matrix(n, n, -11, 24)
            B = generate_random_matrix(n, n, -11, 24)

            expected = np.dot(A, B)
            given, flops = run_ai(A, B)

            print(n, flops)

            # Compare each element of the matrix
            for a in range(A.shape[0]):
                for b in range(A.shape[1]):
                    self.assertEqual(np.round(expected[a][b], 5), np.round(given[a][b], 5))


def generate_random_matrix(rows, columns, min_val, max_val):
    matrix = np.zeros((rows, columns))
    for a in range(rows):
        for b in range(columns):
            matrix[a][b] = randint(min_val, max_val)

    return np.array(matrix)


if __name__ == '__main__':
    unittest.main()
