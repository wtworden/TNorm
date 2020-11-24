
import os
import shutil

from .sage_types import pi, sin, cos, arccos, arcsin, arctan, sqrt, span, vector, Matrix, n

VERTEX_RADIUS = '.03'
AXIS_RADIUS = '.008'
CONE_RADIUS = '.03'
CONE_HEIGHT = '.09'
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

def make_x3d_text(text,coords,scale,size,color):
    text = "\n<Transform translation='{0}' scale='{1}'><Billboard axisOfRotation='0 0 0'><Shape><Appearance><Material ambientIntensity='0.0933' diffuseColor='{4} {5} {6}' shininess='0.51' specularColor='{4} {5} {6}'></Material></Appearance><Text string='{2}' solid='false'><FontStyle family="'Times' 'Orbitron'" quality='3.0' size='{3}'></FontStyle></Text></Shape></Billboard></Transform>".format(coords,scale,text,size,color[0],color[1],color[2])
    return text

def make_x3d_label(genus,bdy,coords,scale,size1,size2,color):
    x=size1*.78
    y=size1*(-.3)
    text = "\n<Transform translation='{0} {1} {2}' scale='{3}'><Billboard axisOfRotation='0 0 0'><Group><Shape><Appearance><Material ambientIntensity='0.0933' diffuseColor='{9} {10} {11}' shininess='0.51' specularColor='{9} {10} {11}'></Material></Appearance><Text string='{6}' solid='false'><FontStyle family="'Times' 'Orbitron'" quality='3.0' size='{4}'></FontStyle></Text></Shape><Transform translation='{12} {13} 0.0'><Shape><Appearance><Material ambientIntensity='0.0933' diffuseColor='{9} {10} {11}' shininess='0.51' specularColor='{9} {10} {11}'></Material></Appearance><Text string='&#160 &#160 {7},{8}' solid='false'><FontStyle family="'Times' 'Orbitron'" quality='3.0' size='{5}'></FontStyle></Text></Shape></Transform></Group></Billboard></Transform>".format(coords[0],coords[1],coords[2],scale,size1,size2,'&#931',genus,bdy,color[0],color[1],color[2],x,y)
    return text

def make_x3d_vertex(coords,radius):
    vert = "\n<Transform translation='{}'><Shape>\n<Appearance>\n<Material ambientIntensity='0.0933' diffuseColor='1 1 1' shininess='0.51' specularColor='0.46 0.46 0.46'></Material>\n</Appearance>\n<Sphere ccw='true' lit='true' radius='{}' solid='true' subdivision='24,24'>\n</Sphere>\n</Shape>\n</Transform>".format(coords,radius)
    return vert

def make_x3d_edge(v1_coords, v2_coords,radius,color):
    v1 = v1_coords
    v2 = v2_coords
    str_v1 = '{} {} {}'.format(rnd(v1[0]),rnd(v1[1]),rnd(v1[2]))
    v = v2-v1
    rho = length = rnd(v.norm())
    trans = length/2
    phi = rnd(arccos(v[2]/rho))
    if abs(v[0]) < 0.0001:
        theta = n(pi/2) if v[1]>0 else n(-pi/2)
    else:
        theta = rnd(arctan(v[1]/v[0]))
        if v[0] < 0:
            theta = rnd(theta + pi)
    edge = "\n<Transform translation='{5}'>\n<Transform rotation='0 0 1 {4}' center='0 0 0'>\n<Transform rotation='0 1 0 {3}' center='0 0 0'>\n<Transform rotation='1 0 0 1.57079' center='0 0 0'>\n<Transform translation='0.0 {2} 0.0'>\n<Shape>\n<Appearance>\n<Material diffuseColor='{6} {7} {8}' shininess='1.0' specularColor='{6} {7} {8}'>\n</Material>\n</Appearance>\n<Cylinder radius='{0}' height='{1}' solid='true'>\n</Cylinder>\n</Shape>\n</Transform>\n</Transform>\n</Transform>\n</Transform>\n</Transform>".format(radius,length,trans,phi,theta,str_v1,color[0],color[1],color[2])
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

#def make_x3d_schlegel(schlegel_projection, scale_factor, color):
#    SP = schlegel_projection
#    edges = SP.lines
#    coords = list(map(lambda v: v*scale_factor,list(map(rnd_vec,SP.coordinates_of(SP.points)))))
#    x3d_str = '\n'
#    x3d_str += "<Shape>\n<Appearance><Material ambientIntensity='0' emissiveColor='{2} {3} {4}' diffuseColor='{2} {3} {4}' specularColor='{2} {3} {4}' shininess='0.0078125' transparency='0'></Material></Appearance><IndexedLineSet coordIndex='{0}' colorPerVertex='false'>\n<Coordinate point='{1}'>\n</Coordinate>\n</IndexedLineSet>\n</Shape>"
#    coordIndex = ''
#    coordPoint = str(coords).replace('0.00000','0.0').replace('[(','').replace(')]','').replace(',','').replace('(','').replace(')',',')
#    for p in SP.lines:
#        coordIndex += '{} {} -1 '.format(p[0],p[1])
#    x3d_str = x3d_str.format(coordIndex, coordPoint, color[0], color[1], color[2])
#    return x3d_str

def make_x3d_schlegel_edges(schlegel_projection, scale_factor, radius, color):
    SP = schlegel_projection
    edges = SP.lines
    coords = list(map(lambda v: v*scale_factor,list(map(rnd_vec,SP.coordinates_of(SP.points)))))
    x3d_str = '\n'
    for line in SP.lines:
        edge = make_x3d_edge(coords[line[0]],coords[line[1]], radius, color)
        x3d_str += edge
    return x3d_str

def make_x3d_html(norm_ball, transparency=.7, online=False, directory=None, dual=False, show_labels=True, show_axes=True, highlight_edges=True, label_color=(1,1,1), edge_color=(1,1,1), label_size=1, edge_thickness=.01, vertex_size=.06, ball_color=(0,0,1), draw_schlegel_faces=False, window_size=(1000,700)):
    if dual:
        title = 'Dual TNorm Ball'
    else:
        title = 'TNorm Ball'
    text = ''
    P=norm_ball.polyhedron()

    width = window_size[0]
    height = window_size[1]

    if P.dim() == 4:

        facets = sorted(P.faces(3), key=lambda x: len(x.vertices()))

        # choose a face with maximal number of vertices. We will project from just outside this face.
        top_facet = facets[-1]
        tf_verts = [top_facet.vertices()[i].vector() for i in range(len(top_facet.vertices()))]
        vecs = [tf_verts[i]-tf_verts[0] for i in range(1, len(tf_verts))]
        N = normal(*span(vecs).basis())
        N_eps = .1*N

        # the center point of our projection face
        center = top_facet.as_polyhedron().center()

        # make sure N is an outward normal, and flip it if its not.
        in_, out_ = P.contains(center+N_eps), P.contains(center-N_eps)
        while not ( in_ or out_ ):
            N_eps *= .5
            in_, out_ = P.contains(center+N_eps), P.contains(center-N_eps)
        if in_:
            N = -N

        # we will rotate P so that N is in the direction of the positive w-axis, then project to the 
        # (x,y,z)-hyperplane. In retrospect this is probably not the most efficient approach. Should 
        # probably skip the rotation, and just project to the plane determined by N and passing through
        # the origin. 
        U = N.normalized()
        phi1, phi2, eta = get_hopf_coords(U,plane='xy')
        phi_xy = pi/2-phi1
        phi1, phi2, eta = get_hopf_coords(vector([0, n(sin(eta)), U[2], U[3]]), plane='yz')
        phi_yz = pi/2-phi1
        phi1, phi2, eta = get_hopf_coords(vector([0, 0, sin(eta), U[3]]), plane='xy')
        phi_zw = pi/2-phi2

        def rotated(vector4D):
            return rotate(vector4D, phi_xy, phi_yz, phi_zw)
        
        vertices = [v.vector() for v in P.vertices()]
        rotated_verts = list(map(rotated,vertices))
        rotated_center = rotated(center)

        # translate the polytope so that the center of the projection face is on the w-axis.
        translation = n(vector([rotated_center[0], rotated_center[1], rotated_center[2],0]))
        translated_verts = list(map(n,[v-translation for v in rotated_verts]))
        translated_center = rotated_center-translation

        # vertex indices of edges, faces, and facets (3d faces)
        edges = [(e.vertices()[0].index(), e.vertices()[1].index()) for e in P.faces(1)]
        faces = [tuple(f.vertices()[i].index() for i in range(len(f.vertices()))) for f in P.faces(2)] 
        facets = [tuple(f.vertices()[i].index() for i in range(len(f.vertices()))) for f in facets] 

        # we need to make sure the vertices of each face are cyclically ordered so x3dom can draw faces 
        fixed_faces = []
        for face in faces:
            face_verts = list(face)
            fixed_face = [face_verts.pop()]
            while len(face_verts) > 0:
                i=0
                while True:
                    if (fixed_face[-1],face_verts[i]) in edges or (face_verts[i],fixed_face[-1]) in edges:
                        fixed_face.append(face_verts[i])
                        face_verts.remove(face_verts[i])
                        break
                    else:
                        i += 1
            fixed_faces.append(fixed_face)
        faces = fixed_faces

        # We project from a point just above the top_facet. This point is chosen to be halfway between
        # the center point of the face (which is now on the w-axis), and the lowest point on the w-axis
        # that is above the center point, and that is in the hyperplane of one of the other facets.
        center_w_int = translated_center[3]
        w_intersections = []
        for f in facets[:-1]:
            vecs = [vertices[f[i]] - vertices[f[0]] for i in range(1,len(f))]
            N = rotated(normal(*span(vecs).basis()))
            k = N.dot_product(translated_verts[f[0]])
            if abs(N[3]) > 0.001:
                w_int = k/N[3]
                if w_int > center_w_int:
                    w_intersections.append(w_int)
        min_w_int = min(w_intersections)
        project_from_w = min((min_w_int + center_w_int)/2.0, center_w_int*2.0)
        project_from = vector([0,0,0,n(project_from_w)])

        # now we project all vertices to the (x,y,z)-hyperplane.
        projected_vertices = []
        for Q in translated_verts:
            PQ = Q - project_from
            t = -project_from[3]/PQ[3]
            projected_vert = project_from+t*PQ
            projected_vertices.append(n(projected_vert[:-1]))
        
        vertex_coords = projected_vertices

        ####### commented out becuase Sage's Schlegel projection does not work correctly as of (11/20/20)
        ####### it looks like this has been fixed and will be working in the next release, so we may want
        ####### to go back to using Sage's built in Schegel projection. So don't delte this...
        #SP = P.schlegel_projection()
        #hull=Polyhedron(SP.coordinates_of(SP.points))
        #coords = list(map(rnd_vec,SP.coordinates_of(SP.points)))
        #bounding_box = (vector([max([v[i] for v in coords]) for i in range(3)]),vector([min([v[i] for v in coords]) for i in range(3)]))

        # get the bounding box and from this a scale factor, which will allow us to put everything at a 
        # scale so default values for the size of edges, vertices, etc will always be reasonable.
        bounding_box = (vector([min([v[i] for v in vertex_coords]) for i in range(3)]),vector([max([v[i] for v in vertex_coords]) for i in range(3)]))
        scale_factor = n(3.0/(max(bounding_box[1])))

        # generate x3dom code for highlighted edges of norm ball
        edges_str = ''
        for e in edges:
            edge = make_x3d_edge(scale_factor*vertex_coords[e[0]],scale_factor*vertex_coords[e[1]], edge_thickness/2., edge_color)
            edges_str += edge
        #edges_str = make_x3d_schlegel_edges(SP, scale_factor, edge_thickness/2., edge_color)
        
        # generate x3dom code for 2D faces of the projected norm ball
        faces_str = ''
        if draw_schlegel_faces:
            coordIndex = ''
            scaled_coords = list(map(lambda v: v*scale_factor,list(map(rnd_vec,vertex_coords))))
            coordPoint = str(scaled_coords).replace('0.00000','0.0').replace('[(','').replace(')]','').replace(',','').replace('(','').replace(')',',')
            for face in faces:
                #vertices = [face.vertices()[i].vector() for i in range(len(face.vertices()))]
                #sp_indices = [SP.coords.index(v) for v in vertices]
                indices_str = ''
                for ind in face:
                    indices_str += '{} '.format(ind)
                indices_str += '-1, '
                coordIndex += indices_str
            coordIndex = coordIndex[:-2] #remove the extra comma at the end

            faces_str = '<Transform scale=\'1 1 1\'>\n<Shape>\n<IndexedFaceSet coordIndex=\'{0}\' solid=\'false\'>\n<Coordinate point=\'{1}\'>\n</Coordinate>\n</IndexedFaceSet>\n<Appearance>\n<Material ambientIntensity=\'.7\' transparency=\'{5}\' diffuseColor=\'{2} {3} {4}\' shininess=\'0.3\' specularColor=\'1.0 1.0 1.0\'>\n</Material>\n</Appearance>\n</Shape>\n</Transform>'.format(coordIndex, coordPoint, ball_color[0], ball_color[1], ball_color[2], transparency)


    else:
        bounding_box = P.bounding_box()
        max_coord = max(bounding_box[1])
        scale_factor = 3./max_coord
        faces_str = P.render_solid().transform(scale=(scale_factor,scale_factor,scale_factor)).x3d_str()
        face_str = fix_x3d_str(faces_str).replace('><','>\n<')
        faces_str = faces_str.replace('Material ','Material ambientIntensity=\'.7\' transparency=\'{}\' '.format(transparency))
        faces_str = faces_str.replace('specularColor=\'0.0 0.0 0.0\'','specularColor=\'1.0 1.0 1.0\'')
        faces_str = faces_str.replace('diffuseColor=\'0.0 0.0 1.0\'','diffuseColor=\'{} {} {}\''.format(ball_color[0], ball_color[1], ball_color[2]))
        faces_str = faces_str.replace('shininess=\'1.0\'','shininess=\'0.3\'')
        vertex_coords = P.vertices_list()

        edges_str = ''
        if highlight_edges:
            for e in P.faces(1):
                v1, v2 = e.vertices()[0].vector(), e.vertices()[1].vector()
                v1, v2 = n(v1,20)*scale_factor, n(v2,20)*scale_factor
                edges_str += make_x3d_edge(v1,v2,edge_thickness/2.,edge_color)

    vertices_str = ''
    labels_str = ''
    for i in range(len(vertex_coords)):
        v = vertex_coords[i]
        scaled_v = n(vector(v),20)*scale_factor
        normalized_v = n(vector(v),20).normalized()
        label_coords = scaled_v + normalized_v*.25
        v = norm_ball.vertex(i)
        if show_labels:
            if dual:
                label = '{}'.format(i)
                str_label_coords = '{} {} {}'.format(label_coords[0], label_coords[1], label_coords[2])
                labels_str += make_x3d_text(label,str_label_coords,'.3 .3 .3',label_size,label_color)
            else:
                labels_str += make_x3d_label(v.genus(),v.num_boundary_comps(),label_coords,'.3 .3 .3',label_size,.7*label_size,label_color)
        str_vert_coords = '{} {} {}'.format(scaled_v[0], scaled_v[1], scaled_v[2])
        vertices_str += make_x3d_vertex(str_vert_coords, VERTEX_RADIUS)


    xmin, xmax = n(bounding_box[0][0]*scale_factor,20), n(bounding_box[1][0]*scale_factor,20)
    ymin, ymax = n(bounding_box[0][1]*scale_factor,20), n(bounding_box[1][1]*scale_factor,20)
    zmin, zmax = n(bounding_box[0][2]*scale_factor,20), n(bounding_box[1][2]*scale_factor,20)

    if show_axes:
        datadir = os.path.dirname(__file__)
        path = os.path.join(datadir, 'data', 'x3d_templates', 'axes.html')
        with open(path,'r') as f:
            slack = min([(xmax-xmin)/2., (ymax-ymin)/2., (zmax-zmin)/2.])
            axis_len = ((xmax-xmin)+slack, (ymax-ymin)+slack, (zmax-zmin)+slack)
            axes = ''.join(f.readlines()).format(axis_len[1], axis_len[1]/2., axis_len[1]/2.+.2, axis_len[0], axis_len[0]/2., axis_len[0]/2.+.2,axis_len[2], axis_len[2]/2., axis_len[2]/2.+.2,AXIS_RADIUS, CONE_RADIUS, CONE_HEIGHT)
    else:
        axes = ''

    if online:
        head = "<head>\n<title>{}</title>\n<script type='text/javascript' src='http://www.x3dom.org/x3dom/release/x3dom.js'>\n</script>\n<link rel='stylesheet' type='text/css' href='http://www.x3dom.org/x3dom/release/x3dom.css'>\n</link>\n<h1>{}</h1>\n<style>\n</style>\n</head>\n".format(title,text)
    else:
        head = "<head>\n<title>{}</title>\n<script type='text/javascript' src=x3dom.js>\n</script>\n<link rel='stylesheet' type='text/css' href=x3dom.css></link>\n<h1>{}</h1>\n<style>\n</style>\n</head>\n".format(title,text) 




#    xaxes = "\n<Transform rotation='1 0 0 0' translation='0 0 0'>\n<Shape>\n<Appearance><Material ambientIntensity='0' emissiveColor='0 0 0' diffuseColor='0 0 0' specularColor='0 0 0' shininess='0.0078125' transparency='0'></Material></Appearance><IndexedLineSet coordIndex='0 1 -1 2 1 -1 3 1 -1' colorPerVertex='false'>\n<Coordinate point='{} 0 0, {} 0 0, '>\n</Coordinate>\n</IndexedLineSet>\n</Shape>\n</Transform>".format(xmin-2,xmax+2,)

    #light = "\n<directionalLight id='directional' direction='0 -1 0' on ='TRUE' intensity='2.0' shadowIntensity='0.0'>\n</directionalLight>"

    body = "<body>\n<x3d width='{}px' height='{}px', margin-left='300px'>\n<scene>\n<Background backUrl='[]' bind='false' bottomUrl='[]' crossOrigin='""' description='""' frontUrl='[]' groundAngle='[]' groundColor='(0.7,0.7,0.75)' isActive='false' leftUrl='[]' rightUrl='[]' scaling='[]' skyAngle='[]' skyColor='(0.7,0.7,0.75)' topUrl='[]' transparency='0/1' ></Background>\n{}\n</scene>\n</x3d>\n</body>\n".format(width,height,faces_str+labels_str+vertices_str+edges_str+axes)

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


def get_hopf_coords(vector, plane='xy'):
    x,y,z,w = n(vector[0]), n(vector[1]), n(vector[2]), n(vector[3])

    Pi = n(pi)

    if plane == 'xy':
        if abs(x) < 0.0001:
            phi = Pi/2 if y>0.000001 else -Pi/2
        else:
            phi = n(arctan(y/x)) + Pi*(x<0)
    
        if phi < 0:
            phi += 2.0*Pi
    
        if abs(z) < 0.0001:
            theta = Pi/2.0 if w>0.000001 else -Pi/2.0
        else:
            theta = n(arctan(w/z)) + Pi*(z<0)
    
        if theta < 0:
            theta += 2.0*Pi
    
        eta = n(arcsin(sqrt(x**2+y**2)))

    elif plane == 'yz':
        if abs(y) < 0.0001:
            phi = Pi/2.0 if z>0.000001 else -Pi/2.0
        else:
            phi = n(arctan(z/y)) + Pi*(y<0)
    
        if phi < 0:
            phi += 2.0*Pi
    
        if abs(x) < 0.0001:
            theta = Pi/2.0 if w>0.000001 else -Pi/2.0
        else:
            theta = n(arctan(w/x)) + Pi*(x<0)
    
        if theta < 0:
            theta += 2.0*Pi
    
        eta = n(arcsin(sqrt(y**2+z**2)))

    return n(phi), n(theta), n(eta)

# rotate vector4D by angle phi_ij about the ij-plane for ij in [xy, zw, xz, yw].
def rotate(vector4D, phi_xy, phi_yz, phi_zw):


    phi1, phi2, eta = get_hopf_coords(vector4D.normalized(), plane='xy')
    r = vector4D.norm()
    rotated_xy = vector([n(cos(phi1+phi_xy)*sin(eta)), n(sin(phi1+phi_xy)*sin(eta)), n(cos(phi2)*cos(eta)), n(sin(phi2)*cos(eta))])
    phi1, phi2, eta = get_hopf_coords(rotated_xy, plane='yz')
    rotated_yz = vector([n(cos(phi2)*cos(eta)), n(cos(phi1+phi_yz)*sin(eta)), n(sin(phi1+phi_yz)*sin(eta)), n(sin(phi2)*cos(eta))])
    phi1, phi2, eta = get_hopf_coords(rotated_yz, plane='xy')
    rotated_zw = vector([n(cos(phi1)*sin(eta)), n(sin(phi1)*sin(eta)), n(cos(phi2+phi_zw)*cos(eta)), n(sin(phi2+phi_zw)*cos(eta))])
    return rotated_zw*r

def translate(vector, translation_vector):
    pass

# get normal vector to the hyperplane containing four given points
#def normal(p1, p2, p3, p4):
#
#    N1 = Matrix([p2-p1, p3-p1, p4-p1, vector([1,0,0,0])]).determinant()
#    N2 = Matrix([p2-p1, p3-p1, p4-p1, vector([0,1,0,0])]).determinant()
#    N3 = Matrix([p2-p1, p3-p1, p4-p1, vector([0,0,1,0])]).determinant()
#    N4 = Matrix([p2-p1, p3-p1, p4-p1, vector([0,0,0,1])]).determinant()
#
#    return vector([N1,N2,N3,N4])
#
# get normal vector to the hyperplane containing four given points
def normal(v1, v2, v3):

    N1 = Matrix([v1, v2, v3, vector([1,0,0,0])]).determinant()
    N2 = Matrix([v1, v2, v3, vector([0,1,0,0])]).determinant()
    N3 = Matrix([v1, v2, v3, vector([0,0,1,0])]).determinant()
    N4 = Matrix([v1, v2, v3, vector([0,0,0,1])]).determinant()

    return vector([N1,N2,N3,N4])



