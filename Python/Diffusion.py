#
# Define Person
#
import numpy as np

class Person(object):
    def __init__( self, pid=0 ):
        self.pid = pid
        self.status = "Potential"

    def __str__( self ):
        msg = '%d,%s' % \
            ( self.pid, self.status )
        return msg

    def adopt(self):
        self.status = "Adopt"

    def dispose(self):
        self.status = "Dispose"


class Population(object):
    """
    simulate a population of Persons
    in a 3-compartment diffusion of innovation model """

    def __init__(self, N=100, beta=0.05, gamma=0.03):
        self.N = N
        self.beta = beta
        self.gamma = gamma

        self.num_potentials = self.N
        self.num_adopters = 0
        self.num_disposers = 0

        self.person_list = []
        self.setup()

    def __str__( self ):
        msg = " "
        msg = '%d,%g,%g,%g,%g,%g' % \
         (self.N,
            self.beta,
            self.gamma,
            self.num_potentials,
            self.num_adopters,
            self.num_disposers)
        return msg

    def setup( self ):
        for i in range( self.N ):
            p = Person( pid=i )
            chance = np.random.random( )
            if chance < self.beta:
                p.adopt( )
                self.num_potentials -= 1
                self.num_adopters += 1
            self.person_list.append( p )

    def model_adoption( self ):
        for i in range( self.N ):
            p = self.person_list[i]
            if ( p.status == "Potential" ):
                chance = np.random.random( )
                if chance < self.beta:
                    p.adopt( )
                    self.num_potentials -= 1
                    self.num_adopters += 1


    def model_disposal( self ):
        for i in range( self.N ):
            p = self.person_list[i]
            if ( p.status == "Adopt" ):
                chance = np.random.random( )
                if chance < self.gamma:
                    p.dispose( )
                    self.num_adopters -= 1
                    self.num_disposers += 1


def main():
    pass

if __name__ == '__main__':
    main()