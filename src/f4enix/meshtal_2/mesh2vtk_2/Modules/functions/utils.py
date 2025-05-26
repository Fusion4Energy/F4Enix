import numpy 

class myOpen():  
    def __init__(self,filename, mode):
        self._file = open(filename, mode)
    
    def readline(self):
        return self._file.readline()
    
    def tell(self):
        return self._file.tell()
   
    def seek(self,pos):
        self._file.seek(pos)

    def skipline(self,n=1,verbose=False):
        for i in range(n):
            line = self._file.readline()
            if verbose : print (f'skipped: {line}')     

    def get_multilines(self,nvalues):
        values = []
        while len(values) < nvalues:
            line = self._file.readline()
            values.extend(line.split())
        return numpy.array(values,dtype=numpy.double)    

class ExtraBin(numpy.ndarray):
    """class storing energy, time or specfic userbin """
    def __new__(cls,data:numpy.ndarray, bintag:str, explicitbin:bool=True ):
        obj = super().__new__(cls,data.shape)
        obj[:] = data[:]
        return obj
    
    def __init__(self, data:numpy.ndarray, bintag:str, explicitbin:bool=True ):
       self._type = bintag
       self._explicit = explicitbin


       if bintag in ('erg','tme'):
            self._binbound = True
            self._totalbin = len(data) != 2
       elif bintag in ('dtme','nuc','cel','line'): 
            self._binbound = False   
            if bintag == 'dtme':
                self._totalbin = False
            else:     
                self._totalbin = len(data) != 1        
       else:  
            raise('bad bin type label')
       
       self._nvalue = len(data)-1 if self._binbound else len(data)

    @property
    def nvalue(self) -> int:
        """
        number of intervals or discrte values 
        """
        return self._nvalue
          
    @property
    def type(self) -> str:
        """
        type of data stored in the bin (erg, tme, nuc, cel, ...)
        """
        return self._type

    @property
    def explicit(self) -> bool:
        """
        True the bin is explcitly defined in the mesh definition.
        """
        return self._explicit

    @property
    def binbound(self) -> bool:
        """
        True if data are bin boundaries.
        """
        return self._binbound

    @property
    def totalbin(self) -> bool:
        """
        True if a total bin is included in data
        """
        return self._totalbin