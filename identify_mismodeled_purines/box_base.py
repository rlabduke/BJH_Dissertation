import sys,os
import math
import numpy as np
import map_value_at_point

def get_dot_product(a,b) :
  assert len(a) == len(b), 'Can\'t take dot of vectors with differing lengths.'
  answer = 0
  for i in range(len(a)) : answer += (a[i]*b[i])
  return answer

def get_determinate_2d(a,b) :
  assert len(a) == len(b)
  assert len(a) == 2
  return (a[0]*b[1]) - (a[1]*b[0])

def get_cross_product(a,b) :
  assert len(a) == len(b)
  assert len(a) == 3
  x = get_determinate_2d((a[1],a[2]),(b[1],b[2]))
  y = -1*(get_determinate_2d((a[0],a[2]),(b[0],b[2])))
  z = get_determinate_2d((a[0],a[1]),(b[0],b[1]))
  return (x,y,z)

def get_amplitude(a) :
  sum = 0
  for i in range(len(a)) : sum += a[i]**2
  return math.sqrt(sum)

def get_unit_vector(v) :
  mag = get_amplitude(v)
  l = [e/mag for e in v]
  return tuple(l)

def subtract_vectors(v1,v2) :
  assert len(v1) == len(v2)
  l = [v1[i] - v2[i] for i in range(len(v1))]
  return tuple(l)

def get_perpindicular_vector_to_plane(p,q,r) :
  assert len(p) == len(q) and len(p) == len(r) and len(p) == 3
  # where p,q,r are points on the plane
  # first find two vectors on the plane, v1 and v2
  # NOTE that the vectors being crossed is pq and pr
  v1 = (p[0] - q[0], p[1] - q[1], p[2] - q[2])
  v2 = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
  return get_cross_product(v1,v2)

def translate_vector(v,new_origin) :
  # translate a vector from the origin to new_origin and return the
  # point of the translated vector, the tain is new_origin
  assert len(v) == len(new_origin)
  return (v[0]+new_origin[0],v[1]+new_origin[1],v[2]+new_origin[2])

def multiply_matricies(m1,m2) :
  # where m1 and m2 is a list(tuples) of lists(tuples)
  assert len(m1[0]) == len(m2), 'column length of m2 does not equal row length of m1'
  m1_rows = m1
  m2_columns = []
  m2_row_len = len(m2[0])
  for col_i in range(m2_row_len) :
    col = [row[col_i] for row in m2]
    m2_columns.append(col)
  answer = []
  for m1_row in m1_rows :
    new_row = []
    for m2_col in m2_columns :
      new_row.append(get_dot_product(m1_row,m2_col))
    answer.append(new_row)
  return tuple(answer)

def scale_vector(v,scale):
  return tuple([e * scale for e in v])

def get_rotation_matrix(axis,theta,angle_format='degree') :
  if angle_format == 'degree' : theta = theta*(math.pi/180)
  elif angle_format != 'radians' :
    raise RuntimeError('angle_format not recognized')
  axis = get_unit_vector(axis)
  a = math.cos(theta/2)
  b,c,d = scale_vector(axis,-1*(math.sin(theta/2)))
  return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                   [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                   [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def rotate_ab_on_abc(a,b,c,theta,angle_format='degree') :
  # rotate ab about the normal to the plane abc centered at a and return head and tail positions
  ba = subtract_vectors(b,a)
  axis = get_perpindicular_vector_to_plane(a,b,c)
  # test perpindicular_vector visually
  # print kin_vector_list([a,translate_vector(axis,a)],'red');sys.exit()
  rm = get_rotation_matrix(axis,theta,angle_format)
  new_vector = get_dot_product(rm,ba)
  # test new_vector visually
  # print kin_vector_list([a,translate_vector(new_vector,a)],'red');sys.exit()
  # print kin_vector_list((ab,(0,0,0),new_vector),'red');sys.exit()
  return a,translate_vector(new_vector,a)

def kin_vector_list(points,color='cyan',group_name='new vector') :
  kin = '@group {%s}\n' % group_name
  kin += '@vectorlist {wqe} color= %s\n{1} P' % color
  for i in range(len(points[0])) :
    kin += ' %.3f' % points[0][i]
  kin += '\n'
  if len(points) == 1 : return kin
  for i in range(len(points) - 1) :
    kin += '{%i} ' % (i + 2)
    for j in range(len(points[i])) :
      kin += ' %.3f' % points[i + 1][j]
    kin += '\n'
  return kin

def get_kin_balls(points,color='cyan',group_name='points') :
  kin = '@group {%s} dominant\n' % group_name
  kin += '@balllist {wtf} color= %s radius = 0.02\n' % color
  for k,v in points.items() :
    id = '%s %s %s' % k
    kin += '{%s} P %.3f %.3f %.3f\n' % (id,v[0],v[1],v[2])
  return kin

def box_base(base,base_type) :
  assert base_type.upper() in ['DA','DT','DC','DG']
  if base_type.upper() in ['DA','DG'] :
    for e in ['C6','N1','C2'] :
      assert e in base.keys()
    dist_from_a2_a3 = 2.5
    side1_scale = 8
    side2_scale = 8.5
    atom_1 = base['C6']
    atom_2 = base['N1']
    atom_3 = base['C2']
  if base_type.upper() in ['DT','DC'] :
    for e in ['N3','C4','C5'] :
      assert e in base.keys()
    dist_from_a2_a3 = 2.5
    side1_scale = 6
    side2_scale = 7
    atom_1 = base['N3']
    atom_2 = base['C4']
    atom_3 = base['C5']
  # FIRST: get perpindicular to base plane
  perp_to_base = get_perpindicular_vector_to_plane(atom_1,atom_2,atom_3)
  # test perp_to_base visually
  # tv = translate_vector(perp_to_base,atom_2)]
  # print kin_vector_list([atom_2,tv,'red');sys.exit()
  # SECOND: get perpindicular to plane defined by atom_2, atom_3 and the 
  # perp_to_base translated origin at atom_2 
  perp_to_b_trans = translate_vector(perp_to_base,atom_2)
  # print kin_vector_list([perp_to_b_trans,atom_2],'green');sys.exit()
  perp_to_perp  = get_perpindicular_vector_to_plane(perp_to_b_trans,atom_2,atom_3)
  # THIRD : scale perp_to_perp, which is essentially perpendicular to the
  # atom_2,atom_3 vector on the base plane
  perp_to_perp_scale = scale_vector(get_unit_vector(perp_to_perp),dist_from_a2_a3)
  # FOURTH : translate the atom_2,atom_3 vector to the end of 
  # perp_to_perp_scale 
  p2p_trans = translate_vector(perp_to_perp_scale,atom_2)
  #   p2p_trans is a point on the base plane perp to atom_2,atom_3
  # FIFTH : get atom_2,atom_3 and do a side1_scale/2 scale to find point_1
  side_1_v = subtract_vectors(atom_2,atom_3)
  side_1_unit = get_unit_vector(side_1_v)
  side_1_v_scaled_half = scale_vector(side_1_unit,side1_scale*(0.5))
  # print kin_vector_list([p2p_trans,atom_2],'green');sys.exit()
  # side_1_v_scaled_half will be the 1/2 of one of the sides of the square.
  # point_2 will be 1/2 the side1_scale up from p2p_trans.
  point_2 = translate_vector(side_1_v_scaled_half,p2p_trans)
  # SIXTH : get the negative of side_1_unit vector scaled to side1_scale to 
  # find second point.
  side_1_v_scaled = scale_vector(side_1_unit,-1*side1_scale)
  point_3 = translate_vector(side_1_v_scaled,point_2)
  # SEVENTH : find the vector perp to side_1_unit and perp_to_base, scale to
  # side2_scale, and translate to point_2 and point_1 get the third and forth
  # points, respectively. 
  side_2_unit = get_unit_vector(get_cross_product(perp_to_base,side_1_unit))
  side_2_scaled = scale_vector(side_2_unit,side2_scale)
  point_4 = translate_vector(side_2_scaled,point_3)
  point_1 = translate_vector(side_2_scaled,point_2)
  # test visually
  # print kin_vector_list([point_1,point_2,point_3,point_4,point_1],'red');sys.exit()
  # returns the corners of the box 1-4 where you have  
  #    1---------2
  #    |         |
  # bb |    b    |
  #    |  /      |
  #    4---------3
  # where b is the base and bb is the backbone
  # print point_1,point_2,point_3,point_4 
  return point_1,point_2,point_3,point_4

def get_scaled_vector_on_line(pA, pB, scale) :
  vec = subtract_vectors(pA,pB)
  vec_unit = get_unit_vector(vec)
  return scale_vector(vec_unit,scale)

class TestPoints(object) :
  def __init__(self, point_1, point_2) :
    self.point_1 = point_1
    self.point_2 = point_2

def get_test_points(base,base_type,sample_spacing) :
  p1,p2,p3,p4 = box_base(base,base_type)
  # point 1
  over_scale = 1.75
  up_scale = 5.5
  up_vec_scaled = get_scaled_vector_on_line(p1, p4, up_scale)
  up_vec = translate_vector(up_vec_scaled, p4)
  over_vec_scaled = get_scaled_vector_on_line(p2, p1, over_scale)
  point_1 = translate_vector(over_vec_scaled, up_vec)
#   print kin_vector_list([p4,up_vec,point_1],'green')
  # point 2
  over_scale = 3.5
  up_scale = 1.25
  up_vec_scaled = get_scaled_vector_on_line(p1, p4, up_scale)
  up_vec = translate_vector(up_vec_scaled, p4)
  over_vec_scaled = get_scaled_vector_on_line(p3, p4, over_scale)
  point_2 = translate_vector(over_vec_scaled, up_vec)
#   print kin_vector_list([p4,up_vec,point_2],'green')
  
  return TestPoints(point_1, point_2)

def get_points_around_base(base,base_type,sample_spacing) :
  p1,p2,p3,p4 = box_base(base,base_type)
  side_1 = subtract_vectors(p2,p1)
  side_2 = subtract_vectors(p4,p1)
  side1_length = get_amplitude(side_1)
  side2_length = get_amplitude(side_2)
  # side_1 goes from p1 to p2
  # side_2 goes from p1 to p4
  #tv1 = translate_vector(side_1,p1)
  #print kin_vector_list([p1,tv1], 'green','s1')
  #tv2 = translate_vector(side_2,p1)
  #print kin_vector_list([p1,tv2], 'green','s2')
  # Now we take the unit vector of side 1 and 2 and scale it to the sample_spacing.
  # We then march down side one, increasing by sample_spacing, taking sample points
  # at scaled side_2 position down.
  side1_unit = get_unit_vector(side_1)
  side2_unit = get_unit_vector(side_2)
  tv1 = translate_vector(side1_unit,p1)
  # print kin_vector_list([p1,tv1], 'green');sys.exit()
  # print side1_length
  # we need to get 2 layers above and below the base plane. For this we need two
  # points above and below p1 perpendicular to the base plane.
  perp_vector = get_unit_vector(get_cross_product(side1_unit,side2_unit))
  perp_vector_scaled = scale_vector 
  p1s = []
  for i in range(-2,3) :
    if i == 0 : p1s.append(p1)
    else :
      #p1s.append(translate_vector(scale_vector(perp_vector,sample_spacing*i),p1))
      p1s.append(translate_vector(scale_vector(perp_vector,0.35*i),p1))
  sample_points = []
  sp = {}
  for pi,p in enumerate(p1s) :
    column = []
    for i in range(1,int(side2_length/sample_spacing) + 1) :
      row = []
      point_on_2 = translate_vector(scale_vector(side2_unit,sample_spacing*i),p)
      for j in range(1,int(side1_length/sample_spacing)) :
        tv = translate_vector(scale_vector(side1_unit,sample_spacing*j),point_on_2)
        row.append(tv)
        sp[(pi,i,j)] = tv
      column.append(row)
    sample_points.append(column)
  # print len(plane), len(plane[0]), len(plane[0][0]);sys.exit()
  sample_points = [e for columns in sample_points for rows in columns for e in rows]
  
  #print get_kin_balls(sp);sys.exit()
  return sp
  return sample_points

def run(args) :
  mtz_file, pdb_file = None, None
  for arg in args :
    if arg.endswith('.pdb') : pdb_file = arg
    elif arg.endswith('.mtz') : mtz_file = arg
  assert pdb_file != None, 'pdb file not provided'
  #assert mtz_file != None, 'mtz file not provided'
  da_base = {}
  da_base['C2'] = (1.424,-58.286,-15.808)
  da_base['N1'] = (1.588,-58.581,-14.509)
  da_base['C6'] = (0.889,-57.872,-13.587)
  points = box_base(da_base,'DA')
  # print points
  #get_points_around_base(da_base,'DA',2.9*(1/4.))
  


  da_base = {}
  da_base['C2'] = (1.424,-58.286,-15.808)
  da_base['N1'] = (1.588,-58.581,-14.509)
  da_base['C6'] = (0.889,-57.872,-13.587)
  # get_points_around_base(da_base,'DA',0.25)
  
  dg_base = {}
  dg_base['c2'] = (11.067,-49.784,-7.982)
  dg_base['n1'] = (10.694,-49.461,-9.266)
  dg_base['c6'] = (9.766,-48.473,-9.595)
  #box_base(dg_base,'G')
  
  dc_base = {}
  dc_base['n3'] = (-1.648,-60.053,-14.445)
  dc_base['c4'] = (-0.8,-61.023,-13.941)
  dc_base['c5'] = (-0.12,-61.834,-14.925)
  # box_base(dc_base,'C')

  dt_base = {}
  dt_base['n3'] = (39.371,-19.689,-10.449)
  dt_base['c4'] = (38.958,-18.852,-9.427)
  dt_base['c5'] = (39.641,-19.035,-8.165)
  #box_base(dt_base,'T')
  

if __name__ == '__main__' :
  run(sys.argv[1:])

