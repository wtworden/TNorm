
import os
import tempfile
import shutil

from .sage_types import *

VERTEX_RADIUS = '.045'
AXIS_RADIUS = '.015'
CONE_RADIUS = '.05'
CONE_HEIGHT = '.15'
EDGE_RADIUS = '.005'

def fix_x3d_str(x3d_str):
    x3d_str = x3d_str.replace('\n','')
    x3d_split = x3d_str.split('>')
    x3d_tags = []
    for string in x3d_split:
        split_str = string.split('<')
        for ss in split_str:
            if ss.replace(' ','') != '':
                x3d_tags.append(ss)
    x3d_str = ''
    for tag in x3d_tags:
        if tag[-1] == '/':
            tag_name = tag.split(' ')[0]
            formatted_tag = '<'+tag[:-1]+'></{}>'.format(tag_name)
        else:
            formatted_tag = '<'+tag+'>'
        x3d_str += formatted_tag
    return x3d_str

def make_x3d_text(text,coords,scale):
    text = "\n<Transform translation='{}' scale='{}'><Billboard axisOfRotation='0 0 0'><Shape><Text string='{}' solid='false'></Text><Appearance><Material diffuseColor='1.0 1.0 1.0' shininess='1.0' specularColor='1.0 1.0 1.0'></Material></Appearance></Shape></Billboard></Transform>".format(coords,scale,text)
    return text

def make_x3d_vertex(coords,radius):
    vert = "\n<Transform translation='{}'><Shape>\n<Sphere ccw='true' lit='true' radius='{}' solid='true' subdivision='24,24'>\n</Sphere>\n<Appearance>\n<Material diffuseColor='1.0 1.0 1.0' shininess='1.0' specularColor='1.0 1.0 1.0'>\n</Material>\n</Appearance>\n</Shape>\n</Transform>".format(coords,radius)
    return vert

def make_x3d_edge(v1_coords, v2_coords,radius):
    v1 = v1_coords
    v2 = v2_coords
    str_v1 = '{} {} {}'.format(rnd(v1[0]),rnd(v1[1]),rnd(v1[2]))
    v = v2-v1
    rho = length = rnd(v.norm())
    trans = length/2
    phi = rnd(arccos(v[2]/rho))
    theta = rnd(arctan(v[1]/v[0]))
    if v[0] < 0:
        theta = rnd(theta + pi)

    #edge = "\n<Transform translation='0.0 {2} 0.0'>\n<Transform rotation='1 0 0 1.57079' center='0 0 0'>\n<Transform rotation='0 1 0 {3}' center='0 0 0'>\n<Transform rotation='0 0 1 {4}' center='0 0 0'>\n<Transform translation {5}><Shape>\n<Cylinder radius='{0}' height='{1}' solid='true'>\n</Sphere>\n<Appearance>\n<Material diffuseColor='1.0 1.0 1.0' shininess='1.0' specularColor='1.0 1.0 1.0'>\n</Material>\n</Appearance>\n</Shape>\n</Transform>\n</Transform>\n</Transform>\n</Transform>\n</Transform>".format(radius,length,trans,phi,theta,str_v1)
    edge = "\n<Transform translation='{5}'>\n<Transform rotation='0 0 1 {4}' center='0 0 0'>\n<Transform rotation='0 1 0 {3}' center='0 0 0'>\n<Transform rotation='1 0 0 1.57079' center='0 0 0'>\n<Transform translation='0.0 {2} 0.0'>\n<Shape>\n<Cylinder radius='{0}' height='{1}' solid='true'>\n</Sphere>\n<Appearance>\n<Material diffuseColor='1.0 1.0 1.0' shininess='1.0' specularColor='1.0 1.0 1.0'>\n</Material>\n</Appearance>\n</Shape>\n</Transform>\n</Transform>\n</Transform>\n</Transform>\n</Transform>".format(radius,length,trans,phi,theta,str_v1)
    return edge    

def rnd(x):
    if abs(x) < .00001:
        return 0
    else:
        return n(x,20)   

def rnd_vec(v):
    return vector(map(rnd,v))

#def scale_vec(v, factor):
#    return v*factor

def make_x3d_schlegel(schlegel_projection, scale_factor):
    SP = schlegel_projection
    edges = SP.lines
    coords = map(lambda v: v*scale_factor,map(rnd_vec,SP.coordinates_of(SP.points)))
    x3d_str = '\n'
    x3d_str += "<Shape>\n<Appearance><Material ambientIntensity='0' emissiveColor='1 1 1' diffuseColor='1 1 1' specularColor='1 1 1' shininess='0.0078125' transparency='0'></Material></Appearance><IndexedLineSet coordIndex='{}' colorPerVertex='false'>\n<Coordinate point='{}'>\n</Coordinate>\n</IndexedLineSet>\n</Shape>"
    coordIndex = ''
    coordPoint = str(coords).replace('0.00000','0.0').replace('[(','').replace(')]','').replace(',','').replace('(','').replace(')',',')
    for p in SP.lines:
        coordIndex += '{} {} -1, '.format(p[0],p[1])
    x3d_str = x3d_str.format(coordIndex, coordPoint)
    return x3d_str


def make_x3d_html(norm_ball, online=False, directory=None, dual=False):
    if dual:
        title = 'Dual TNorm Ball'
    else:
        title = 'TNorm Ball'
    text = ''
    P=norm_ball.polyhedron()

    if P.dim() == 4:
        SP = P.schlegel_projection()
        hull=Polyhedron(SP.coordinates_of(SP.points))
        bounding_box = hull.bounding_box()
        scale_factor = 3/(max(bounding_box[1]))
        x3d_str = make_x3d_schlegel(SP, scale_factor)
        vertices_coords = map(lambda x: n(x,20),SP.coordinates_of(SP.points))

    else:
        bounding_box = P.bounding_box()
        max_coord = max(bounding_box[1])
        scale_factor = 3./max_coord
        x3d_str = P.render_solid().transform(scale=(scale_factor,scale_factor,scale_factor)).x3d_str()
        x3d_str = fix_x3d_str(x3d_str).replace('><','>\n<')
        x3d_str = x3d_str.replace('Material ','Material ambientIntensity=\'.7\' transparency=\'.4\' ')
        x3d_str = x3d_str.replace('specularColor=\'0.0 0.0 0.0\'','specularColor=\'1.0 1.0 1.0\'')
        x3d_str = x3d_str.replace('shininess=\'1.0\'','shininess=\'0.3\'')
        vertices_coords = P.vertices_list()

    vertices_str = ''
    labels_str = ''
    for i in range(len(vertices_coords)):
        v = vertices_coords[i]
        scaled_v = n(vector(v),20)*scale_factor
        normalized_v = n(vector(v),20).normalized()
        label_coords = scaled_v + normalized_v*.25
        str_vert_coords = '{} {} {}'.format(scaled_v[0], scaled_v[1], scaled_v[2])
        str_label_coords = '{} {} {}'.format(label_coords[0], label_coords[1], label_coords[2])
        v = norm_ball.vertex(i)
        if dual:
            label = '{}'.format(i)
        else:
            frac = '1/{}*'.format(-v.euler_char()) if v.euler_char()!=-1 else ''
            label = '{}S_{},{}'.format(frac,v.genus(),v.num_boundary_comps())
        vertices_str += make_x3d_vertex(str_vert_coords, VERTEX_RADIUS)
        labels_str += make_x3d_text(label,str_label_coords,'.3 .3 .3')
    edges_str = ''
    for e in P.faces(1):
        v1, v2 = e.vertices()[0].vector(), e.vertices()[1].vector()
        v1, v2 = n(v1,20)*scale_factor, n(v2,20)*scale_factor
        edges_str += make_x3d_edge(v1,v2,EDGE_RADIUS)

    xmin, xmax = n(bounding_box[0][0]*scale_factor,20), n(bounding_box[1][0]*scale_factor,20)
    ymin, ymax = n(bounding_box[0][1]*scale_factor,20), n(bounding_box[1][1]*scale_factor,20)
    zmin, zmax = n(bounding_box[0][2]*scale_factor,20), n(bounding_box[1][2]*scale_factor,20)

    datadir = os.path.dirname(__file__)
    path = os.path.join(datadir, 'data', 'x3d_templates', 'axes.html')
    with open(path,'r') as f:
        axis_len = max(xmax-xmin+2, ymax-ymin+2, zmax-zmin+2)
        axes = ''.join(f.readlines()).format(axis_len, axis_len/2., axis_len/2.+.2, AXIS_RADIUS, CONE_RADIUS, CONE_HEIGHT)


    if online:
        head = "<head>\n<title>{}</title>\n<script type='text/javascript' src='http://www.x3dom.org/x3dom/release/x3dom.js'>\n</script>\n<link rel='stylesheet' type='text/css' href='http://www.x3dom.org/x3dom/release/x3dom.css'>\n</link>\n<h1>{}</h1>\n<style>\n</style>\n</head>\n".format(title,text)
    else:
        head = "<head>\n<title>{}</title>\n<script type='text/javascript' src=x3dom.js>\n</script>\n<link rel='stylesheet' type='text/css' href=x3dom.css></link>\n<h1>{}</h1>\n<style>\n</style>\n</head>\n".format(title,text) 




#    xaxes = "\n<Transform rotation='1 0 0 0' translation='0 0 0'>\n<Shape>\n<Appearance><Material ambientIntensity='0' emissiveColor='0 0 0' diffuseColor='0 0 0' specularColor='0 0 0' shininess='0.0078125' transparency='0'></Material></Appearance><IndexedLineSet coordIndex='0 1 -1 2 1 -1 3 1 -1' colorPerVertex='false'>\n<Coordinate point='{} 0 0, {} 0 0, '>\n</Coordinate>\n</IndexedLineSet>\n</Shape>\n</Transform>".format(xmin-2,xmax+2,)

    #light = "\n<directionalLight id='directional' direction='0 -1 0' on ='TRUE' intensity='2.0' shadowIntensity='0.0'>\n</directionalLight>"

    body = "<body>\n<x3d width='1000px' height='700px', margin-left='300px'>\n<scene>\n<Background backUrl='[]' bind='false' bottomUrl='[]' crossOrigin='""' description='""' frontUrl='[]' groundAngle='[]' groundColor='(0.7,0.7,0.75)' isActive='false' leftUrl='[]' rightUrl='[]' scaling='[]' skyAngle='[]' skyColor='(0.7,0.7,0.75)' topUrl='[]' transparency='0/1' ></Background>\n{}\n</scene>\n</x3d>\n</body>\n".format(x3d_str+labels_str+vertices_str+edges_str+axes)

    page = "<html>\n{}{}</html>".format(head,body)

    datadir = os.path.dirname(__file__)
    path = os.path.join(datadir, 'data', 'x3dom')
    css = os.path.join(path, 'x3dom.css')
    js = os.path.join(path, 'x3dom.js')

    shutil.copy2(css, directory)
    shutil.copy2(js, directory)

    filename = os.path.join(directory, 'normball_x3d.html')
    
    with open(filename,'w') as file:
        file.write(page)
    return os.path.abspath(filename)















