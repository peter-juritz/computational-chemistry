class Basis:
    def read_basis(self,fname):
        f = open(fname,'r')
        state = 0 # 0: waiting for atom 1: waiting for func  2: reading func
        bs = [ ]
        k=0
        n = 0
        ats ={}
        cur = ''
        while (True):
            str = f.readline()
            sp = str.split()
            if (str == ''):
                break
            if (sp == []):
                continue
            if(sp[0]=='****'):
                bs = []
                state = 0
                continue
            if (sp[0][0] == '!'):
                continue
            if (state == 0):
                cur = sp[0].lower()
                ats[cur] = []
                ats[cur].append(sp[0])
                state=1
                continue
            elif (state == 1):
                if (sp[0] == 'S'):
                    state = 2
                elif (sp[0] == 'SP'):
                    state = 3
                ats[cur].append([sp[0]])
                n = int(sp[1])
                continue
            elif (state == 2):
                ats[cur][len(ats[cur])-1].append([float(sp[0]),float(sp[1])]) # TODO improve
                n-=1
                if (n==0):
                    state = 1
                    continue
            elif(state == 3):
                ats[cur][len(ats[cur])-1].append([float(sp[0]),float(sp[1]),float(sp[2])])      # XXX hack
                n-=1
                if (n==0):
                    state = 1
        return ats



