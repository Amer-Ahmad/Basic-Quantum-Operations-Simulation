import random
import cmath

# this if the start of the first class for the Quantum register,
# the class basically is a fancy vector, so every object instantiated out of this class is a vector
# additionally with a set of operations that you can execute on this vector (methods).
class QuantumRegister:

    # constructor, takes the default value of N as 1 (N represents the number of qubits in the system), if the amplitudes array was not passed as a parameter,
    # then it will set this array to ket{0...0} which is quivalent to the vector[1,0,0,...0]
    # (self is like this in java)
    def __init__(self, N=1, amplitudes=None ):
        self.N = N
        if amplitudes is None:
            # setting the amplitudes array if it was not passed as a parameter
            amplitudes = [complex(0, 0) for i in range(2 ** N)]
            amplitudes[0] = complex( 1 , 0 )
        self.amplitudes = amplitudes


    # this method is to calculate the norm of any Quantum register, although the norm will always be equal to 1, implementing
    # this method is for the sake of validating Quantum states that we have.
    def norm(self):
        ret = 0
        # the norm gets calculated by squaring the probability amplitudes (the element in the vector of the quantum state)
        # in complex numbers r^2 = x^2 + y^2.
        for x in self.amplitudes:
            ret += ( x.imag ** 2 + x.real ** 2 )
        return ret

    # this method is to calculate the inner product between two Quantum registers, which will be equal to scalar.
    # we calculate the inner product by multiplying the first vector with the conjugunt of the second vector.
    # note: the conjugunt of a complex number z = x + yi is equal to z` = x - yi
    def InnerProduct(self, Q2):
        ret = 0
        for x,y in zip(self.amplitudes, Q2.amplitudes):
            conj = complex( x.real , -1*x.imag )
            ret += ( conj * y )
        return ret

    # this method is to calculate the tensor product between two Quantum registers, and it produces a quantum register
    # with dimensionality of n1+n2 where n1 is the dimension of Q1, and n2 is the dimension of Q2,
    # we multiply the two vectros like we are multiplying polynomials.
    def TensorProduct(self, Q2):
        # Q3 is the resulting Quantum register.
        Q3 = QuantumRegister( self.N + Q2.N )
        for x in range (len(self.amplitudes)):
            for y in range (len(Q2.amplitudes)):
                Q3.amplitudes[ y + x*len(Q2.amplitudes) ] = self.amplitudes[x] * Q2.amplitudes[y]
        return Q3

    # this method does a measurement to the Quantum register based on its probability amplitudes,
    # and then it returns the equivalent computational basis in decimal (e.g. if we got |101> after the measurment we
    # the method will return 5)
    def measurement(self):
        # generating a random number so we do the measurment.
        # we previously know that the probability of choosing any of the computational basis is equal to its probability
        # amplitude squared.
        # note: after a measurement we have to get one of the computational basis.
        rand = random.random()
        sum=0
        state=-1
        decimalEquivalent=-1
        #in this loop we are just finding which state we are got after the measurement.
        for x in range (len(self.amplitudes)):
            sum += ( self.amplitudes[x].imag**2 + self.amplitudes[x].real**2 )
            if rand<sum:
                ans = self.amplitudes
                # saving the decimal quivalent of the quantum state that i got after the measurement.
                decimalEquivalent = x
                break

        # setting up the quantum register of the state that i got after the measurment.
        for x in range (2 ** self.N):
            if x == decimalEquivalent:
                self.amplitudes[x] = 1
            else:
                self.amplitudes[x] = 0

        return decimalEquivalent


    # a method that takes a permutation vector and permute the current quantum state accoring to the permutation vector.
    # it's self explanatory.
    def permute(self, permutation):
        temp = [None] * len(permutation)
        for x in range (len(permutation)):
            temp[x] = self.amplitudes[ permutation[x]-1 ]
        self.amplitudes = temp

# this class is like a fancy matrix, additionally with operations to apply on it.
class LinearOperator:

    # constructor, N default value is equal to 1, in case the matrix of the linear operator was not passed as a parameter.
    # then we set the matrix to the identity matrix.
    def __init__(self, N=1, matrix=None ):
        self.N = N

        # in case of not passing the matrix as a parameter.
        if matrix is None:
            matrix = [None] * ( 2 ** N )
            for x in range (2 ** N):
                matrix[x] = [None] * (2 ** N)

            for x in range (2 ** N):
                for y in range (2 ** N):
                    if x is y:
                        matrix[x][y] = complex(1.0 , 0)
                    else:
                        matrix[x][y] = complex(0.0 , 0)
        self.matrix = matrix

    # this is another constructor that initialize the linear operator matrix as the outer product of two Quantum states.
    # as the out of product of a vector Nx1 with another vetor Mx1 gives the matrix NxM
    @classmethod
    def initialize(cls, f, g):
        tmp = [None] * (2 ** f.N)
        for x in range(2 ** f.N):
            tmp[x] = [None] * (2 ** f.N)

        for x in range(2 ** f.N):
            for y in range(2 ** g.N):
                #print( f.amplitudes[x] * g.amplitudes[y] )
                tmp[x][y] = f.amplitudes[x] * g.amplitudes[y]

        return cls(f.N, tmp)

    # this method add two linear operators with each other, in addition of having the option to multiply each of the linear
    # operators with a scalar.
    def add(self, L2, alpha, beta):

        ret = [None] * (2 ** self.N)
        tmp = [None] * (2 ** self.N)
        for x in range(2 ** self.N):
            ret[x] = [None] * (2 ** self.N)
            tmp[x] = [None] * (2 ** self.N)

        for x in range(2 ** self.N):
            for y in range(2 ** self.N):
                ret[x][y] = self.matrix[x][y] * alpha
                tmp[x][y] = L2.matrix[x][y] * beta
                ret[x][y] += tmp[x][y]

        return ret


    # this method takes a Quantum register u, and then it multiplies u with the linear operator,
    # basically applying a function on the quantum state u.
    def transform(self, u):
        tmp = [0] * (2 ** self.N)
        for x in range(2 ** self.N):
            for y in range(2 ** self.N):
                tmp[x] += self.matrix[x][y] * u.amplitudes[y]
        u.amplitudes = tmp

    # this method takes a Quantum state that has a dimensionality bigger than the linear operator, and it takes also a vector
    # that contains certain indices to apply the transformation on.
    # so basically the same as the previous method except on certain indices.
    def partial_transform(self, u, indices):
        tmp = [0] * (2 ** self.N)
        for x in range(2 ** self.N):
            for y in range(2 ** self.N):
                tmp[x] += self.matrix[x][y] * u.amplitudes[ indices[y] ]
        u.amplitudes = tmp


# for testing I'm initializing some amplitude arrays taking into consideration that the norm of the vector of any quantum state
# will always be equal to 1.
amplitude = [ complex(0.5 , 0), complex(1/cmath.sqrt(2) , 0), complex(1/cmath.sqrt(2) , 0), complex(0.0 , 0) ]
amplitude4 = [ complex(0.25 , 0), complex(0.25 , 0), complex(0.5 , 0), complex( cmath.sqrt(5)/cmath.sqrt(8) , 0) ]
amplitude5 = [ complex(0.25 , 0), complex(0.25 , 0), complex(0.0 , 0), complex( 0.25 , 0) , complex( 0.0 , 0) , complex( 0.0 , 0 ), complex(0.25 , 0) , complex( 0.0 , 0 ) ]
Q1 = QuantumRegister(2 , amplitude)
Q4 = QuantumRegister(2 , amplitude4 )
Q2 = QuantumRegister(2)
Q5 = QuantumRegister(3 , amplitude5 )
# print( Q1.norm() )
# print( Q1.InnerProduct(Q2) )
# Q3 = Q1.TensorProduct(Q2)
# print(Q3.amplitudes)
#
# print( Q1.measurement() )
# print( Q1.amplitudes )
#
# print( Q4.norm() )
# Q4.permute( (2 , 1 , 4 , 3) )
# print( Q4.amplitudes )

#testing the linear operators.

L1 = LinearOperator(2)
print(L1.matrix)

L2 = LinearOperator.initialize(Q1, Q4)
print(L2.matrix)

print( L1.add( L2, complex(1.0 , 0) , complex( 2.0 , 0 ) ) )

L2.transform(Q4)
print(Q4.amplitudes)

L2.partial_transform(Q5, (0, 1, 3, 6))
print( Q5.amplitudes )