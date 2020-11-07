import matplotlib.pyplot as plt
import Diffusion as df

class DOI_Model(object):
    def __init__( self, N=500, beta=0.09, gamma=0.01, max_time=250 ):

        self.potential_history = []
        self.adoption_history = []
        self.disposal_history = []

        self.N = N
        self.beta = beta
        self.gamma = gamma
        self.max_time = max_time
        self.timespan = range(max_time)

        self.pop = df.Population(N, beta, gamma)

    def __str__( self ):
        msg = 'DOI_Model: %d, %8.2f, %8.2f, %d' % \
              (self.N, self.beta, self.gamma, self.max_time)
        return msg

    def run(self):
        for time in self.timespan:
            # keep track of what happened
            self.potential_history.append(self.pop.num_potentials)
            self.adoption_history.append(self.pop.num_adopters)
            self.disposal_history.append(self.pop.num_disposers)

            # model population actions
            self.pop.model_adoption()
            self.pop.model_disposal()

    def plot(self):

        # plot results
        #
        plt.plot(self.timespan, self.potential_history, '-b', label="potential")
        plt.plot(self.timespan, self.adoption_history, '-r', label="adopters")
        plt.plot(self.timespan, self.disposal_history, '-g', label="disposers")

        plt.title('diffusion of innovation with random mixing')
        plt.xlabel('time')
        plt.legend(title="key", loc='center right')
        plt.ylabel('adoption rate')

        plt.ylim(0, pop.N)

        plt.show()


#
# object-oriented main #
def main():
    model = DOI_Model(N=500, beta=0.08, gamma=0.02, max_time = 250)
    model.run()
    model.plot()
#
# use standard Python idiom #
if __name__ == '__main__':
    main()