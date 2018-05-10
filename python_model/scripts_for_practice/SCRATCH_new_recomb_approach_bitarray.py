import bitarray
import numpy as np
import numpy.random as r


#TODO:
    #rename stuff (paths??? and this isn't a queue! call it what it is)
    #try out integrating it into the model (maybe its own branch), test-running it, and timing it



#define a class that will:
    #a. hold a large, low-mem numpy array of recomb paths, to be used as a queue of recomb paths
    #b. hold a backup bitarray of recomb paths, which
        #the numpy array will be refreshed from every time it runs out
    #c. hold a pointer to the recomb_paths bitarray file, which 
        #the backup bitarray will be reloaded from every time it runs out
class Recomb_Queue:
    def __init__(self, r_lookup, queue_len, bitarray_len, queue_len_min = 1000, recomb_file = './.recomb.tmp', bitarray_shuffle_multiple = 100, n_recomb_per_loop = 10000):
        self.r_lookup = r_lookup
        self.queue_len = queue_len 
        self.queue_len_min = queue_len_min
        self.recomb_file = recomb_file
        self.L = r_lookup.shape[0]
        self.bitarray_len = bitarray_len
            #NOTE: bitarray will write 0s to end of file if bitstring not divisible by 8, so this calculates the amount of bits to remove when file is read back in
        self.drop_file_bits = -1*(self.bitarray_len*self.L)%8
        self.bitarray_len_min = self.L * self.queue_len_min
        self.bitarray_shuffle_multiple = bitarray_shuffle_multiple
        self.n_recomb_per_loop = n_recomb_per_loop
        
        self.bitarray = create_recomb_path_bitarray(self.r_lookup, n_paths = bitarray_len, recomb_file = self.recomb_file, n_recomb_per_loop = n_recomb_per_loop)
        
        self.queue = None
        self.refresh_queue(start = True) 
        

    def refresh_queue(self, start = False):
        print('\nREFRESHING QUEUE\n')
        tot = self.L * self.queue_len
        new_queue = np.reshape(np.int8(self.bitarray[:tot]), (self.L, self.queue_len), order = 'F')
        if start == False:
            self.queue = np.hstack((self.queue, new_queue))
            self.bitarray = self.bitarray[tot:]
            if len(self.bitarray) < self.bitarray_len_min:
                #read in the bitarray file again 
                self.refresh_bitarray()
        else:
            self.queue = new_queue
            self.bitarray = self.bitarray[tot:]


    def refresh_bitarray(self):
        print('\nREFRESHING BITARRAY\n')
        new = bitarray.bitarray()
        with open(self.recomb_file, 'rb') as f:
            new.fromfile(f)
        if self.drop_file_bits != 0:
            new = new[:self.drop_file_bits]
        new = shuffle_bitarray(new, tract_len = self.L * self.bitarray_shuffle_multiple)
        self.bitarray = self.bitarray + new


    def update_queue(self, n):
        self.queue = self.queue[:,n:]

    #get the stipulated number of recombination paths out of the queue
    def get_paths(self, n):
        paths = self.queue[:, :n]
        self.queue = self.queue[:, n:]
        if self.queue.shape[1] < self.queue_len_min:
            self.refresh_queue()
        return(paths)


#generate an array of recombination paths, using the r_lookup object
def gen_recomb_paths(r_lookup, n_recombinants):
    paths = np.array([r.choice(r_lookup[i,], size = n_recombinants, replace = True) for i in range(len(r_lookup))])
    paths = np.cumsum(paths, axis = 0)%2
    paths = list(paths.T.flatten())
    return(paths)



#save the bitarray of recombination paths to a tmp file in the cwd
def create_recomb_path_bitfile(paths_bitarray, recomb_file):
    paths.tofile(recomb_file)


#generate a bitarray of recombination paths, using the gen_recomb_paths function
def create_recomb_path_bitarray(r_lookup, n_paths, recomb_file, n_recomb_per_loop = 10000):
    n_paths = (n_paths//n_recomb_per_loop + int(n_paths%n_recomb_per_loop > 0))*n_recomb_per_loop
    paths = bitarray.bitarray([])
    for it in range(n_paths//n_recomb_per_loop):
        print('\nCreating_array_number %i' % it)
        paths += bitarray.bitarray(gen_recomb_paths(r_lookup, n_recombinants = n_recomb_per_loop))
    create_recomb_path_bitfile(paths, recomb_file)
    return(paths)





#shuffle the genome-length tracts of a bitarray
def shuffle_bitarray(ba, tract_len):
    print('\nSHUFFLING BITARRAY\n')
    assert len(ba)%tract_len == 0, 'ERROR: tract_len is not a factor of len(ba)'
    chunks  = list(range(len(ba)//tract_len))
    r.shuffle(chunks)
    new = bitarray.bitarray([])
    for chunk in chunks: 
        new += ba[chunk*tract_len: (chunk+1)*tract_len]
    return(new)


