import math

def gamma_function(alpha):
    if alpha == 1:
        return 1
    else:
        return (alpha - 1) * gamma_function(alpha - 1)

def gamma_integral(alpha):
    def integrand(t):
        return math.exp(-t) * (t ** (alpha - 1))

    # Approximating the integral using Simpson's rule
    n = 1000  # Number of intervals
    a = 0     # Lower limit of integration
    b = 100   # Upper limit of integration
    h = (b - a) / n

    integral = (integrand(a) + integrand(b)) / 3

    for i in range(1, n):
        integral += (2 + 2 * (i % 2)) * integrand(a + i * h) / 3

    return integral * h

def F(z, m):
    K = (m + 1) / (math.sqrt(m) * math.factorial(m))
    integral = gamma_integral((m + 1) / 2)
    return K * integral

def main():
    degrees_of_freedom = int(input("Enter degrees of freedom (m): "))
    z_values = [float(input(f"Enter z value {i + 1}: ")) for i in range(3)]

    for z in z_values:
        result = F(z, degrees_of_freedom)
        print(f"For z={z}, probability is {result:.6f}")

if __name__ == "__main__":
    main()
