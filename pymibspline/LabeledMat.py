
import sys
import numpy as np
import copy
import time

class LabeledMat:
    """
    Object to store matrix-formatted genetic / genomic data, indexed by row / column names
    
    Usage:
        >>> x = LabeledMat.loadFile(textMatrix.csv, sep=',', dt=a4)

    the `dt` argument stands for data type. `a4` means 4-char strings. You can set them to 
    float, int as well.

    """
    def __init__(self, X, rownames, colnames, verbose=False, log=sys.stderr):
        self.data = X
        self.__len__ = len(self.data)
        self.nrow = X.shape[0]
        self.ncol = X.shape[1]
        self.rownames = rownames
        self.colnames = colnames
        self.rowmap = {}
        self.colmap = {}
        for i in range(len(rownames)):
            self.rowmap[rownames[i]] = i
        for j in range(len(colnames)):
            self.colmap[colnames[j]] = j
        self.log = log
        self.verbose=verbose
        self.infoHeader = "[" + self.__class__.__name__ + "]"

    def __getitem__(self, index):
        if isinstance(index[0], str):
            x = self.rowmap[index[0]]
        elif (isinstance(index[0], list) or isinstance(index[0], np.ndarray)) and isinstance(index[0][0], str):
            x = map(lambda a: self.rowmap[a], index[0])
        elif isinstance(index[0], slice):
            x = [ii for ii in xrange(*index[0].indices(self.nrow))]
        else:
            x = index[0]
            if not isinstance(x, list) and not isinstance(x, np.ndarray) and x < 0:
                x = range(self.nrow)[x]

        if isinstance(index[1], str):
            y = self.colmap[index[1]]
        elif (isinstance(index[1], list) or isinstance(index[1], np.ndarray)) and isinstance(index[1][0],str):
            y = map(lambda a: self.colmap[a], index[1])
        elif isinstance(index[1], slice):
            y = [ii for ii in xrange(*index[1].indices(self.ncol))]
        else:
            y = index[1]
            if not isinstance(y, list) and not isinstance(y, np.ndarray) and y < 0:
                y = range(self.ncol)[y]

        if type(x) == int:
            x = [x]
        if type(y) == int:
            y = [y]
        
        if self.verbose: 
            self.info("idx_x = \n" + str(x))
            self.info("idx_y = \n" + str(y))
        
        return LabeledMat(
                X = self.data[x,][:,y].reshape(len(x), len(y)), 
                rownames = map(lambda a: self.rownames[a], x), 
                colnames = map(lambda a: self.colnames[a], y),
                verbose=self.verbose,
                log = self.log
                )

    def __setitem__(self, index, value):
        if type(index[0]) == str:
            x = self.rowmap[index[0]]
        else:
            x = index[0]
        if type(index[1]) == str:
            y = self.colmap[index[1]]
        else:
            y = index[1]
        self.data[x,y] = value

    def info(self, message):
        print >> self.log, time.strftime("%Y-%m-%d %X ") + self.infoHeader + message

    def transpose(self):
        self.data = self.data.transpose()
        self.nrow = self.data.shape[0]
        self.ncol = self.data.shape[1]
        tmp = self.rownames
        self.rownames = self.colnames
        self.colnames = tmp
        tmp = self.rowmap
        self.rowmap = self.colmap
        self.colmap = tmp

    @staticmethod
    def rbind(X, Y, verbose=False, log=sys.stderr):
        data = np.concatenate((X.data, Y.data), axis=0)
        rownames = copy.deepcopy(X.rownames)
        colnames = copy.deepcopy(X.colnames)
        r2 = copy.deepcopy(Y.rownames)
        
        conflict = set(rownames) & set(r2)
        if len(conflict)>0:
            for c in conflict:
                r2[r2.index[c]] += "_"

        rownames.extend(r2)
        return LabeledMat(data, rownames, colnames, verbose=verbose, log=log)

    @staticmethod
    def cbind(X, Y, verbose=False, log=sys.stderr):
        data = np.concatenate((X.data, Y.data), axis=1)
        rownames = copy.deepcopy(X.rownames)
        colnames = copy.deepcopy(X.colnames)
        c2 = copy.deepcopy(Y.colnames)
        
        conflict = set(colnames) & set(c2)
        if len(conflict)>0:
            for c in conflict:
                c2[c2.index[c]] += "_"
        
        colnames.extend(c2)
        return LabeledMat(data, rownames, colnames, verbose=verbose, log=log)

    @staticmethod
    def loadFile(filename, sep='\t', dt=float, verbose=False, log=sys.stderr):
        fo = open(filename, 'r')
        line = fo.readline()
        tokens = line.strip().split(sep)
        colnames = tokens
        rownames = []
        X = []
        ncols = len(colnames)
        nrows = 0
        firstLine = True
        # first pass: get rownames
        for line in fo:
            tokens = line.strip().split(sep)
            if firstLine:
                firstLine = False
                if len(tokens) == ncols:
                    ncols -= 1
                    colnames = colnames[1:]
            else:
                assert( (len(tokens) - 1 ) == ncols)
            rownames.append(tokens[0])
            nrows += 1
        if verbose:
            print >> log, "Number of rows: \t%s\nNumber of columns: \t%s" % (nrows, ncols)
        fo.close()
        X = np.ndarray(shape=(nrows, ncols), dtype=np.dtype(dt) )
        fo = open(filename, 'r')
        fo.readline()
        rcnt = 0
        if verbose:
            print >> log, "Loading data ..." 
        for line in fo:
            tokens = line.strip().split(sep)
            tokens.pop(0)
            X[rcnt,:] = tokens
            rcnt += 1
        fo.close()
        if verbose:
            print >> log,  "Done."
        return LabeledMat(X, rownames, colnames, verbose=verbose, log=log)

    def write2file(self, filename, sep='\t', eol='\n'):
        fo = open(filename, 'w')
        fo.write("#HEADER" + sep + sep.join(self.colnames) + eol)
        rcnt = 0
        for l in self.data:
            fo.write(self.rownames[rcnt] + sep + sep.join([str(x) for x in l]) + eol)
            rcnt += 1
        fo.close()
        


        
