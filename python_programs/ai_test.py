import unittest
import numpy as np
from random import randint
from ai_recursive import run_ai


class AITest(unittest.TestCase):
    def test_MultiplicationTest(self):
        for n in range(2, 3):
            A = np.array([
                [0.5928, 0.8443, 0.8579, 0.8473, 0.6236, 0.3844, 0.2975, 0.0567, 0.2727, 0.4777],
                [0.8122, 0.4800, 0.3928, 0.8361, 0.3374, 0.6482, 0.3682, 0.9572, 0.1404, 0.8701],
                [0.4736, 0.8009, 0.5205, 0.6789, 0.7206, 0.5820, 0.5374, 0.7586, 0.1059, 0.4736],
                [0.1863, 0.7369, 0.2165, 0.1352, 0.3241, 0.1497, 0.2223, 0.3865, 0.9026, 0.4500],
                [0.6131, 0.9023, 0.0993, 0.9698, 0.6531, 0.1709, 0.3582, 0.7507, 0.6078, 0.3250],
                [0.0384, 0.6343, 0.9589, 0.6528, 0.6351, 0.9953, 0.5818, 0.4144, 0.4747, 0.6235],
                [0.3380, 0.6748, 0.3172, 0.7783, 0.9496, 0.6625, 0.0136, 0.6228, 0.6737, 0.9719],
                [0.8782, 0.5096, 0.0557, 0.4512, 0.0200, 0.4417, 0.9796, 0.3594, 0.4809, 0.6887],
                [0.8805, 0.9182, 0.2168, 0.5652, 0.8651, 0.5090, 0.9167, 0.9212, 0.0831, 0.2777],
                [0.0094, 0.8423, 0.6472, 0.8414, 0.2647, 0.3978, 0.5528, 0.1649, 0.3698, 0.1464]
            ])

            B = np.array([
                [0.9972, 0.9326, 0.1281, 0.9990, 0.2361, 0.3966, 0.3879, 0.6697, 0.9355, 0.8463],
                [0.3133, 0.5245, 0.4435, 0.2296, 0.5344, 0.9140, 0.4572, 0.4307, 0.9391, 0.7784],
                [0.7160, 0.8028, 0.0928, 0.5182, 0.8650, 0.8291, 0.8296, 0.2731, 0.0592, 0.6705],
                [0.5931, 0.6717, 0.4118, 0.1976, 0.2896, 0.1421, 0.7833, 0.4125, 0.0342, 0.6240],
                [0.6606, 0.2985, 0.4461, 0.2221, 0.0734, 0.4692, 0.0962, 0.9034, 0.1195, 0.5248],
                [0.0836, 0.9169, 0.9104, 0.2989, 0.5844, 0.5659, 0.6139, 0.9565, 0.2610, 0.2310],
                [0.5334, 0.9499, 0.4931, 0.5406, 0.7655, 0.0453, 0.1400, 0.7924, 0.0298, 0.8831],
                [0.5408, 0.4480, 0.8921, 0.3776, 0.5384, 0.6523, 0.3613, 0.5710, 0.6378, 0.1263],
                [0.6902, 0.6477, 0.3539, 0.7632, 0.3565, 0.7528, 0.8813, 0.0117, 0.4981, 0.0738],
                [0.7870, 0.0641, 0.3553, 0.9418, 0.3798, 0.7629, 0.7716, 0.3014, 0.7727, 0.1529]
            ])

            expected = np.dot(A, B)
            given, operations_all, mul = run_ai(A, B)

            print(n)

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


def print_matrix(matrix, name="AI_result"):
    print(f"{name}:")

    for row in matrix:
        print("    ".join(f"{elem:.4f}" for elem in row))


if __name__ == '__main__':
    unittest.main()