import numpy as np
import math

# from https://stackoverflow.com/questions/39512260/calculating-gini-coefficient-in-python-numpy
# Gini is the mean absolute difference of all value pairs, normalized by the value average
def gini_coefficient(x):
    """Compute Gini coefficient of array of values"""
    diffsum = 0
    for i, xi in enumerate(x[:-1], 1):
        diffsum += np.sum(np.abs(xi - x[i:]))
    return diffsum / (len(x)**2 * np.mean(x))

def main():
    print("Testing gini_coefficient")
    assert gini_coefficient(np.array([5, 5, 5, 5])) == 0
    # from https://goodcalculators.com/gini-coefficient-calculator/
    assert math.isclose(gini_coefficient(np.array([1000, 2000, 3000, 4000, 5000, 6000, 7000])), 0.28571, rel_tol=0.0001)
    print("gini_coefficient tests passed")


if __name__ == "__main__":
    main()