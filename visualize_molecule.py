from molecular_geometry import GeomReader
from visual import sphere, curve, color
import sys

SIZE_SCALE = 0.005
at_colours = {'c':color.green,'o':color.red,'n':color.blue,'h':color.white,'s':color.black,'fe':color.red,'p':color.blue,'he':color.blue,'f':color.yellow}
vw_radius = {'c':170,'h':120,'o':152,'n':155,'s':180,'fe':240,'p':180,'he':140,'f':147}
def draw_atoms(mol):
    for at in mol.atoms:
    #    if at.name =='h':
            #continue
        atom_sphere = sphere(pos = at.xyz,radius = SIZE_SCALE*(vw_radius[at.name]),color=at_colours[at.name])
        at.vis = atom_sphere

def draw_atoms_vispy(mol):
    from vispy import scene
    from vispy.visuals.transforms import STTransform
    canvas = scene.SceneCanvas(keys='interactive', bgcolor='black',
                            size=(800, 600), show=True)
    view = canvas.central_widget.add_view()
    view.camera = 'arcball'
    for at in mol.atoms:
        atom_sphere = scene.visuals.Sphere(radius=SIZE_SCALE*(vw_radius[at.name]), method='ico', parent=view.scene,
                                    color=at_colours[at.name], edge_color='black',subdivisions=2)

        atom_sphere.transform = STTransform(translate=at.xyz)
    canvas.app.run()

def draw_bonds(mol):
    for l in mol.bonds:
        b = curve(pos =(mol.atoms[l[0]-1].xyz[0],mol.atoms[l[0]-1].xyz[1],mol.atoms[l[0]-1].xyz[2]))
        b.append(mol.atoms[l[1]-1].xyz)

def draw_spatial_center(mol):
    xyz = [0,0,0]
    tot = 0
    for at in mol.atoms:
        tot+= at.mass

    for at in mol.atoms:
        c = at.mass
        xyz[0] += at.x*c
        xyz[1] += at.y*c
        xyz[2] += at.z*c
    xyz[0]/= tot
    xyz[1]/= tot
    xyz[2]/= tot
    sphere(pos= xyz,radius=0.5,color=color.blue)


if __name__ == '__main__':
    if len(sys.argv) not in [2, 3]:
        print 'Usage: %s [method] (geometry_file)' % sys.argv[0]
        print ' Optionally specify a rendering method. Options: vpython, vispy. vpython is default'
        sys.exit(0)
    if len(sys.argv) == 2:
        geom = GeomReader(sys.argv[1])
        draw_atoms(geom.mol)
    if len(sys.argv) == 3:
        geom = GeomReader(sys.argv[2])
        if sys.argv[1] == 'vpython':
            draw_atoms(geom.mol)
        elif sys.argv[1] == 'vispy':
            draw_atoms_vispy(geom.mol)
        else:
            print 'Unsupported method: %s' % sys.argv[1]

