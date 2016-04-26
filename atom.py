from motion import Motion
class Atom:

    name = 'undef';
    def __init__(self,name,number,x=0,y=0,z=0,connect=-1):
        eneg = {'c':2.55,'h':2.1,'o':3.44,'n':3.04}
        molecular_mass = {'c':12.01,'h':0.008,'n':14.01,'o':16}
        self.x = x;
        self.y = y;
        self.z = z;
        self.xyz = [x,y,z];
        self.connect = [connect] # TODO make ptr to mol parent and atoms
        self.name = name
        self.atom_number = number
        self.aof_functions = []
        self.vis = None
        self.v=[0,0,0]
        self.a=[0,0,0]
        masses = {'h':1.0,'he':2.0,'c':6.0,'o':8.0,'n':7,'f':9,'s':16,'p':15}
        self.mass =masses[name]


