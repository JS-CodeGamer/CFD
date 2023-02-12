#! /bin/env python


# Packaged Solver
class solver_packages:
    def backward_difference(diff:int=None, initial_values:list[list[int]]=None):
        return Solver(initial_values, solver_config(-1,diff))

    def central_difference(diff:int=None, initial_values:list[list[int]]=None):
        return Solver(initial_values, solver_config(0,diff))

    def forward_difference(diff:int=None, initial_values:list[list[int]]=None):
        return Solver(initial_values, solver_config(1,diff))


# methods:
#  -1 for backward difference
#   0 for central difference
#   1 for forward difference
# diff: mesh cell size
class solver_config:
    def __init__(self, method:int=None, diff:int=None, layer_size:tuple[int]=None):
        self.method = method
        self.diff = diff
        self.layer_size = layer_size



class Solver:

    def __init__(self, initial_values:list[list[int]]=None, config:solver_config=None):
        if initial_values==None:
            raise TypeError("No Initial conditions given.")
        self.initial_values=initial_values
        self.differentials=[initial_values]
        if config==None or config.diff==None:
            raise TypeError("Solver Configiruation not defined correctly mesh cell size is needed.")
        self.config=config
        self.__load_config()


    def make_new_layer(self):
        new_layer = [ self.config.layer_size[0]*[-1] ]*self.config.layer_size[1]
        self.differentials.append(new_layer)


    def __load_config(self):
        if self.config.method==None:
            print("Warning: Method not defined for solver defaulting to central difference 1st order.")
            self.method=self.__c_diff
        if self.config.method == -1:
            self.method=self.__b_diff
        elif self.config.method == 0:
            self.method=self.__c_diff
        elif self.config.method == 1:
            self.method=self.__f_diff
        if self.config.layer_size==None:
            self.config.layer_size=(len(self.differentials[-1][-1]), len(self.differentials[-1]))


    def get_nth_diff(self, n:int, x:int, y:int):

        # Handle if memory has less layers
        while len(self.differentials)<n :
            self.make_new_layer()

        return self.method(x,y)


    def __b_diff(self, n:int, x:int, y:int):
        # Apply backward difference if this meshpoint has not been visited
        if self.differentials[n][x][y]==-1 :
            # Handle enpoints of mesh (x-coordinate)
            if x==0 :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x+1,y)-self.__b_diff(n-1,x,y))/self.config.diff
            else:
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x-1,y))/self.config.diff

            # Handle enpoints of mesh (y-coordinate)
            if y==0 :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y+1)-self.__b_diff(n-1,x,y))/self.config.diff
            else:
                self.differentials[n][x][y] += (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x,y-1))/self.config.diff

        return self.differentials[n][x][y]


    def __c_diff(self, n:int, x:int, y:int):
        # Apply central difference if this meshpoint has not been visited
        if self.differentials[n][x][y]==-1 :
            # Handle enpoints of mesh (x-coordinate)
            if x==0 :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x+1,y)-self.__b_diff(n-1,x,y))/self.config.diff
            elif x==len(self.differentials[-1]) :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x-1,y))/self.config.diff
            else:
                self.differentials[n][x][y] = (self.__c_diff(n-1,x+1,y)-self.__c_diff(n-1,x-1,y))/(2*self.config.diff)

            # Handle enpoints of mesh (y-coordinate)
            if y==0 :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y+1)-self.__b_diff(n-1,x,y))/self.config.diff
            elif y==len(self.differentials[-1][-1]) :
                self.differentials[n][x][y] += (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x,y-1))/self.config.diff
            else:
                self.differentials[n][x][y] += (self.__c_diff(n-1,x,y+1)-self.__c_diff(n-1,x,y-1))/(2*self.config.diff)

        return self.differentials[n][x][y]


    def __f_diff(self, n:int, x:int, y:int):
        # Apply forward difference if this meshpoint has not been visited
        if self.differentials[n][x][y]==-1 :
            # Handle enpoints of mesh (x-coordinate)
            if x==len(self.differentials[-1]) :
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x-1,y))/self.config.diff
            else:
                self.differentials[n][x][y] = (self.__f_diff(n-1,x+1,y)-self.__f_diff(n-1,x,y))/self.config.diff

            # Handle enpoints of mesh (y-coordinate)
            if y==len(self.differentials[-1][-1]) :
                self.differentials[n][x][y] += (self.__b_diff(n-1,x,y)-self.__b_diff(n-1,x,y-1))/self.config.diff
            else:
                self.differentials[n][x][y] = (self.__b_diff(n-1,x,y+1)-self.__b_diff(n-1,x,y))/self.config.diff

        return self.differentials[n][x][y]


