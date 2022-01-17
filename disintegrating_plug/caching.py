from joblib import Memory
import cloudpickle
from pathlib import Path
import os.path

location = './cachedir'
memory = Memory(location, verbose=2)

class Luggage:

    def __init__(self, fpath):
        self.fpath = fpath
        Path(fpath).mkdir(parents=True, exist_ok=True)

    def memory(self, func):
        fname = self.fpath+'/'+func.__name__+'.pkl'
        def hijack():
            if not os.path.isfile(fname):
                cloudpickle.dump(func(),
                                 open(fname,'wb'))
            return cloudpickle.load(open(fname,'rb'))
        return hijack
