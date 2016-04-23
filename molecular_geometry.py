from  molecule import Molecule
from atom import Atom
import math

from util import (my_dot, scal_dot, vec_add, vec_minus, vec_cross, vec_dot, vec_times, dist, normalise)
class GeomReader:

    def read_cartesian(self,fname):
        molc=Molecule()
        f = open(fname,'r')
        while (True):
            str = f.readline()
            if (str==''):
                break
            sp = str.split()
            if (len(sp) != 4):
                continue
            molc.add_atom(Atom(sp[0].lower(),0,float(sp[1]),float(sp[2]),float(sp[3])))
        self.mol = molc
        return  molc

    def read_zmatrix(self,fname):
        f = open(fname,'r')
        self.file = f;
        name = f.readline().strip()
        self.mol = Molecule()
        self.mol.add_atom(Atom(name.lower(),1,0,0,0,connect=0));

        str = f.readline();
        sp = str.split();
        if (sp == []):
            return self.mol
        sp[0] = sp[0].lower()
        if (sp==[]):
            return self.mol
        try:
            bl = float(sp[2])
        except:
            bl = self.lookup_var(sp[2])
        self.mol.add_atom(Atom(sp[0],2,0,0,bl,int(sp[1])))
        self.mol.add_bond(2,int(sp[1]))



        str = f.readline();
        sp = str.split();
        if (sp == []):
            return self.mol
        sp[0] = sp[0].lower()
        try:
            bl = float(sp[2]);
        except:
            bl = self.lookup_var(sp[2])

        try:
            ba = float(sp[4]);
        except:
            ba = self.lookup_var(sp[4])
        a = Atom(sp[0],3,math.sin(ba*math.pi/180.0)*bl,0,self.mol.atoms[int(sp[1]) -1].z - math.cos(ba*math.pi/180.0)*bl,int(sp[1]))# BIG CHANGE

        self.mol.add_atom(a)
        self.mol.add_bond(3,int(sp[1]))

        atom_number = 4;
        while (True):
            #print atom_number
            str = f.readline();
            sp = str.split();
            if (sp==[]):
                print 'end of matrix'
                break
            sp[0] = sp[0].lower()
            try:
                bl = float(sp[2]);
            except:
                bl = self.lookup_var(sp[2])

            try:
                ba = float(sp[4]);
            except:
                ba = self.lookup_var(sp[4])
            try:
                dh = float(sp[6]);
            except:
                dh = self.lookup_var(sp[6])

            connect = int(sp[1])
            angle_connect = int(sp[3])
            dihed_connect = int(sp[5])


            atoms = self.mol.atoms;
            #print 'connect: %i angle_connect: %i' %(connect,angle_connect)
            vector1 = vec_minus(atoms[connect-1].xyz,atoms[angle_connect-1].xyz );
            vector2 = vec_minus(atoms[connect-1].xyz,atoms[dihed_connect-1].xyz );
            norm1 = vec_cross(vector1,vector2)
            norm2 = vec_cross(vector1,norm1)
            norm1 = normalise(norm1)
            norm2 = normalise(norm2)

            norm1 =vec_times(norm1,-1*math.sin(dh*math.pi/180))
            norm2 = vec_times(norm2,math.cos(dh*math.pi/180))

            vector3 =vec_add(norm1,norm2)
            vector3 =normalise(vector3)

            vector3 = vec_times(vector3,bl*math.sin(ba*math.pi/180.0))

            vector1 = normalise(vector1)

            vector1 = vec_times(vector1,bl*math.cos(ba*math.pi/180.0))

            vector2 = vec_add(atoms[connect - 1].xyz,vector3)
            vector2 = vec_minus(vector2,vector1)

            a = Atom(sp[0],atom_number,vector2[0],vector2[1],vector2[2],int(sp[1]))
            self.mol.add_atom(a)
            self.mol.add_bond(atom_number,int(sp[1]))

            atom_number+=1;
        return self.mol


    def lookup_var(self,name):
       # print 'Looking up ',name
        mark = self.file.tell();
        name.strip()
        mul = 1;
        if (name[0] == '-'):
            mul = -1;
            name = name[1:]

        while (True):
            str = self.file.readline();
            if (str == ''):
                print 'Lookup error'
                break;
            sp = str.split()
            if (sp ==[]):
                continue
            if (sp[0]==name):
                self.file.seek(mark);
                return float(sp[1])*mul;
        print 'Couldnt look up ', name;
        self.file.seek(mark);
        return None;

    def __init__(self, fname):
        if fname.endswith(".zmat"):
            self.read_zmatrix(fname)
        else:
            self.read_cartesian(fname)

    def write_atoms(self,fname):
        f = open(fname,'w')
        for at in self.mol.atoms:
            f.write('%s %f %f %f\n'%(at.name,at.x,at.y,at.z))
        f.close()


    def print_mol_ats(self):
        key = {}
        for at in self.mol.atoms:
            try:
                key[at.name]+=1
            except:
                key[at.name]=1
        for i in key:
            print '%s %i' % (i,key[i])


