import numpy as np
import functools

class DemographicModel:
    '''Stores piecewise-exponential demographic models.'''

    def __init__(self, filename=None):
        # Number of epochs. Must equal the lengths of the following lists:
        self.num_epochs = 0
        # Epoch start times
        self.times = []
        # Epoch starting (i.e. most recent) sizes
        self.sizes = []
        # Epoch growth rates (forward in time)
        self.rates = []

        # If a fastNeutrino output file is specified, read it.
        if filename is not None:
            self.read_fastNeutrino_output(filename)

    def read_fastNeutrino_output(self, model_fn):
        '''Read epochs from a fastNeutrino fitted parameters output file.'''
        with open(model_fn) as modelfile:
            header = modelfile.readline()
            n_anc = float(modelfile.readline())
            # First epoch implicitly starts at t=0
            start_time = 0.0
            for line in modelfile:
                # Get epoch parameters.
                if line.startswith('c'):
                    # Constant-N epoch
                    n, t = map(float, line.split()[-2:])
                    g = 0.0
                elif line.startswith('e'):
                    # Exponential-growth epoch
                    n, t, g = map(float, line.split()[-3:])
                else:
                    raise ValueError('Warning, bad line: ' + line.strip())
                    break

                # Add epoch to model.
                # TODO: handle time order errors.
                self.add_epoch(start_time, n, g)
                # Set next epoch start time to current epoch end time
                start_time = t
        # Add ancestral population size as the last epoch
        self.add_epoch(start_time, n_anc, 0.0)


    def add_epoch(self, time, size, rate):
        '''Add new epoch to the demographic model.'''
        if self.num_epochs == 0 or time > self.times[-1]:
            self.num_epochs += 1
            self.times.append(time)
            self.sizes.append(size)
            self.rates.append(rate)
        else:
            raise ValueError('New epoch time is less than previous epoch time.')

    def population_size(self, T):
        '''Return the population size at time T'''

        if self.num_epochs == 0:
            # No epochs have been added.
            return np.empty_like(T)
        else:
            # Evalute piecewise exponential growth during each epoch.
            conditions = [T>=t0 for t0 in self.times]
            functions = [functools.partial(exponential_growth, *args)
                         for args in zip(self.sizes, self.times, self.rates)]
            return np.piecewise(T, conditions, functions)

def exponential_growth(n0, t0, r, T):
    '''Calculate N(t) for exponentially-growing population back in time'''
    return n0*np.exp(-(T-t0)*r)
