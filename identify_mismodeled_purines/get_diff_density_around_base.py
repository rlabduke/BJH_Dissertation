# (jEdit options) :folding=explicit:collapseFolds=1:
# 
#  This script takes a pdb and mtz to identify difference densities around a 
#  DNA base.
import sys,os
from phenix.automation.refinement import refinement_base
#from mmtbx.command_line import maps
import iotbx.pdb
import mmtbx.utils
from iotbx import crystal_symmetry_from_any
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import mmtbx.maps
from cctbx import crystal
from libtbx import group_args
import libtbx.easy_run
import box_base
import map_value_at_point
import numpy


# {{{ usage
def usage() :
  s = '''
*******************************************************************************
phenix.python get_diff_density_around_base.py xxxx [otiions]

Where xxxx is a pdb id.

To use options type the name of the option then and = sign and the option value
with NO SPACES: i.e. positive_threshold=2.5

Options:
    positive_threshold        The positive sigma value where difference density
                              points at or greater than that value will be 
                              flagged as positive peaks. Default = 2.5
    negative_threshold        The negative sigma value where difference density
                              points at or less than that value will be 
                              flagged as negative peaks. Default = -2.5
    threshold                 The sigma value, both negative and positive, where
                              difference density points <= the negative of that 
                              value or >= that value will be flagged as negative
                              and positive peaks, respectively. Default = 2.5
    point_n_threshold         The smallest number of points to be considered a 
                              difference density peak.
    write_sample_points_kin   Set to True to write a kinemage (xxxx_sample_points.kin)
                              in the current directory that has the sample 
                              points where the fo-fc density was sampled.
*******************************************************************************
'''
  sys.stderr.write(s)
  sys.exit()

# }}}

# {{{ default_params
default_params = """\
maps {
  map_coefficients {
    map_type = 2mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = 2FOFCWT
    mtz_label_phases = PH2FOFCWT
    fill_missing_f_obs = False
  }
  map_coefficients {
    map_type = 2mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = 2FOFCWT_fill
    mtz_label_phases = PH2FOFCWT_fill
    fill_missing_f_obs = True
  }
  map_coefficients {
    map_type = mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = FOFCWT
    mtz_label_phases = PHFOFCWT
    fill_missing_f_obs = False
  }
  map_coefficients {
    map_type = anomalous
    format = *mtz phs
    mtz_label_amplitudes = ANOM
    mtz_label_phases = PHANOM
  }
  map {
    map_type = 2mFo-DFc
    fill_missing_f_obs = False
    grid_resolution_factor = 1/4.
  }
}
"""
# }}}

# {{{ GetMapCoeffs
class GetMapCoeffs(refinement_base):

  def __init__(self, args=None, out=sys.stderr) :
    super(GetMapCoeffs, self).__init__(args=args, out=out)

  def run(self):
    self.get_inputs()
    self.extract_model()
    self.extract_data()

  # {{{ compute_map
  def compute_map(self) :
    log = sys.stderr
    # the mtz and pdb files are self.mtz_file self.pdb_file, respectively
    #out = maps.run(args = [self.pdb_file, self.mtz_file])
    #self.map_coeff_file, self.xplor_maps = out[0],out[1]
    master_params = mmtbx.maps.maps_including_IO_master_params()
    master_params = master_params.fetch(iotbx.phil.parse(default_params))
    processed_args = mmtbx.utils.process_command_line_args(
      args=[self.mtz_file, self.pdb_file],
      log=log,
      master_params=master_params)
    working_phil = processed_args.params
    params = working_phil.extract()
    if params.maps.input.pdb_file_name is None : 
      params.maps.input.pdb_file_name = self.pdb_file
    if params.maps.input.reflection_data.file_name is None :
      params.maps.input.reflection_data.file_name = self.mtz_file
    pdb_inp = iotbx.pdb.input(file_name = params.maps.input.pdb_file_name)
    # get all crystal symmetries
    cs_from_coordinate_files = [pdb_inp.crystal_symmetry_from_cryst1()]
    csa =crystal_symmetry_from_any.extract_from(params.maps.input.reflection_data.file_name)
    cs_from_reflection_files = [csa]
    crystal_symmetry = None
    try :
      crystal_symmetry = crystal.select_crystal_symmetry(
        from_coordinate_files=cs_from_coordinate_files,
        from_reflection_files=cs_from_reflection_files)
    except AssertionError, e :
      if ("No unit cell and symmetry information supplied" in str(e)) :
        raise Sorry("Missing or incomplete symmetry information.  This program "+
          "will only work with reflection file formats that contain both "+
          "unit cell and space group records, such as MTZ files.")
    #
    rfn = params.maps.input.reflection_data.file_name
    rf = reflection_file_reader.any_reflection_file(
        file_name = rfn)
    reflection_files = [rf]
    reflection_file_names = [rfn]
    reflection_file_server = reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry   = True,
      reflection_files = reflection_files, #[],
      err              = log)
    #
    reflection_data_master_params = mmtbx.utils.data_and_flags_master_params(
      master_scope_name="reflection_data")
    reflection_data_input_params = processed_args.params.get(
      "maps.input.reflection_data")
    reflection_data_params = reflection_data_master_params.fetch(
      reflection_data_input_params).extract().reflection_data
    #
    determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
      reflection_file_server  = reflection_file_server,
      parameters              = reflection_data_params,
      data_parameter_scope    = "maps.input.reflection_data",
      flags_parameter_scope   = "maps.input.reflection_data.r_free_flags",
      data_description        = "Reflection data",
      keep_going              = True,
      log                     = log)
    f_obs = determine_data_and_flags_result.f_obs
    self.resolution_range = f_obs.resolution_range()
    r_free_flags = determine_data_and_flags_result.r_free_flags
    test_flag_value = determine_data_and_flags_result.test_flag_value
    if(r_free_flags is None):
      r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
      test_flag_value=None
    print >> log, "-"*79
    print >> log, "\nInput model file:", params.maps.input.pdb_file_name
    pdb_hierarchy = pdb_inp.construct_hierarchy(set_atom_i_seq=True)
    atom_selection_manager = pdb_hierarchy.atom_selection_cache()
    # get xray_structure
    xray_structure = pdb_inp.xray_structure_simple(
      crystal_symmetry = crystal_symmetry)
    mmtbx.utils.setup_scattering_dictionaries(
      scattering_table = params.maps.scattering_table,
      xray_structure   = xray_structure,
      d_min            = f_obs.d_min(),
      log              = log)
    xray_structure.show_summary(f = log, prefix="  ")
    print >> log, "-"*79
    print >> log, "Bulk solvent correction and anisotropic scaling:"
    fmodel = mmtbx.utils.fmodel_simple(
      #update_f_part1_for      = "map",
      xray_structures         = [xray_structure],
      scattering_table        = params.maps.scattering_table,
      f_obs                   = f_obs,
      r_free_flags            = r_free_flags,
      outliers_rejection      = params.maps.input.reflection_data.outliers_rejection,
      skip_twin_detection     = params.maps.skip_twin_detection,
      bulk_solvent_correction = params.maps.bulk_solvent_correction,
      anisotropic_scaling     = params.maps.anisotropic_scaling)
    fmodel_info = fmodel.info()
    fmodel_info.show_rfactors_targets_scales_overall(out = log)
    print >> log, "-"*79
    print >> log, "Compute maps."
    cmo = mmtbx.maps.compute_map_coefficients(
      fmodel = fmodel,
      params = params.maps.map_coefficients,
      pdb_hierarchy = pdb_hierarchy,
      log = log)
    file_name_base = params.maps.input.pdb_file_name
    if(file_name_base.count(".")>0):
      file_name_base = file_name_base[:file_name_base.index(".")]
    self.map_coeff_file_name = file_name_base+"_map_coeffs.mtz"
    r_free_flags_output = None
    write_mtz_file_result = cmo.write_mtz_file(file_name = self.map_coeff_file_name,
      r_free_flags=r_free_flags_output)
    print >> log, "All done."
    if (write_mtz_file_result) :
      print >> log, "Map coefficients: %s" % os.path.basename(self.map_coeff_file_name)
  # }}}

# }}}

# {{{ DensitySamplePoints
class DensitySamplePoints(object) :
  
  def __init__(self) : 
    self.sample_points = {}
  
  def deposit_sample_point(self, grp_arg):
    # grp_arg elements : plane, column, row, xyz, map_value
    k = (grp_arg.plane, grp_arg.column, grp_arg.row)
    self.sample_points[k] = grp_arg

  # {{{ identify_diff_peaks
  def identify_diff_peaks(self, 
                          positive_threshold, 
                          negative_threshold, 
                          count_threshold) :
    # This will identify difference peaks in the difference density around 
    # each base. A peak is defind as at least count_threshold adjacent points of 
    # positive_threshold and above or negative_threshold and below.

    # peak points
    negative_peak_points = {}
    positive_peak_points = {}
    for k,point in self.sample_points.items() :
      if point.map_value <= negative_threshold : 
        negative_peak_points[k] = point
      elif point.map_value >= positive_threshold : 
        positive_peak_points[k] = point

    # identify peaks as defined by function parameters
    sample_matrix = [(i,j,k) for i in range(-1,2) 
                             for j in range(-1,2) 
                             for k in range(-1,2)]
    # sample_matrix.remove((0,0,0))
    self.negative_peaks = []
    for k,point in negative_peak_points.items() :
      # see if point is already in a peak
      peak_or_false = self.point_in_peaks_list(point,self.negative_peaks)
      if peak_or_false == False : 
        peak = {k:point}
        # self.negative_peaks.append(peak)
      else : peak = peak_or_false
      # get adjacent negative_peak_points
      peak_exists = False
      for t in sample_matrix:
        nk = (k[0]+t[0], k[1]+t[1], k[2]+t[2])
        if nk in negative_peak_points.keys() :
          np = negative_peak_points[nk]
          if np not in peak : peak[nk] = np
          npeak_or_false = self.point_in_peaks_list(np,self.negative_peaks)
          if npeak_or_false != False : 
            peak_exists = True
            existingpeak = npeak_or_false
      if peak_exists : 
        for k,v in peak.items() : existingpeak[k] = v
      elif peak_or_false == False : self.negative_peaks.append(peak)

    self.positive_peaks = []
    for k,point in positive_peak_points.items() :
      # see if point is already in a peak
      peak_or_false = self.point_in_peaks_list(point,self.positive_peaks)
      if peak_or_false == False : 
        peak = {k:point}
        # self.positive_peaks.append(peak)
      else : peak = peak_or_false
      # get adjacent positive_peak_points
      peak_exists = False
      for t in sample_matrix:
        nk = (k[0]+t[0], k[1]+t[1], k[2]+t[2])
        if nk in positive_peak_points.keys() :
          np = positive_peak_points[nk]
          if np not in peak : peak[nk] = np
          npeak_or_false = self.point_in_peaks_list(np,self.positive_peaks)
          if npeak_or_false != False : 
            peak_exists = True
            existingpeak = npeak_or_false
      if peak_exists : 
        for k,v in peak.items() : existingpeak[k] = v
      elif peak_or_false == False : self.positive_peaks.append(peak)
    
    # Eliminate peaks with < count_threshold
    for peak in self.negative_peaks :
      if len(peak) < count_threshold : self.negative_peaks.remove(peak)
    for peak in self.positive_peaks :
      if len(peak) < count_threshold : self.positive_peaks.remove(peak)
      
  # }}}

  # {{{ point_in_peaks_list
  def point_in_peaks_list(self, point, peak_list) :
    for peak in peak_list :
      if point in peak.values() : return peak
    return False
  # }}}

# }}}

# {{{ Peak
class Peak(object) :

  # {{{ __init__
  def __init__(self,peak,base) :
    # peak is a grp_args with plane, column, row, map_value and xyz
    self.peak = peak
    self.base = base
    self.n_sample_points = len(peak)
  # }}}

  # {{{ get_3d_distance
  def get_3d_distance(self, x, y) :
    assert len(x) == 3 and len(y) == 3
    sm = 0
    for i in range(3) : sm += (x[i]-y[i])**2
    return numpy.sqrt(sm)
  # }}}

  # {{{ adjacent_to
  def adjacent_to(self,xyz_key,radius=1.5) :
    # returns true if at least one sample ponit
    # is within 'radius' of the given 'xyz_key'
    err = 'xyz_key %s given not found'
    if xyz_key.strip().startswith('test_') : xyz_key = xyz_key.strip()
    else : xyz_key = xyz_key.strip().upper()
    assert xyz_key in vars(self.base.xyz), err % xyz_key
    for k,point in self.peak.items() :
      d = self.get_3d_distance(vars(self.base.xyz)[xyz_key], point.xyz)
      if d <= radius :
        return True
    return False
  # }}}

# }}}

# {{{ BaseClass
class BaseClass(object) :

  # {{{ __init__
  def __init__(self, chain_id, residue, altloc, resolution) :
    self.chain_id   = chain_id
    self.residue    = residue
    self.altloc     = altloc
    self.resolution = resolution
    self.density_sample_points = DensitySamplePoints()
    self.set_atoms_xyz()
    self.set_sample_points()
  # }}}

  # {{{ is_purine, is_pyrimidine
  def is_purine(self) :
    if self.residue.resname.strip().upper() in ['DA','DG'] : return True
    return False

  def is_pyrimidine(self) :
    if self.residue.resname.strip().upper() in ['DC','DT'] : return True
    return False
  # }}}

  # {{{ get_string_id
  def get_string_id(self,delimiter=',') :
    d = delimiter
    s = '%s%s%s%s%s%s%s%s%s'
    return s % (self.chain_id,d,
                self.residue.resseq,d,
                self.residue.icode,d,
                self.altloc,d,
                self.residue.resname)
  # }}}

  # {{{ set_atom_xyz
  def set_atoms_xyz(self, atom_names=None) :
    # atom_names should be a list of string atom names, if None gets all atoms
    # returns a group_args with atom_name->(x,y,z)
    self.xyz = group_args()
    for atom in self.residue.atoms():
      if atom_names == None :
        vars(self.xyz)[atom.name.strip().upper()] = atom.xyz
      elif atom.name.strip() in atom_names :
        vars(self.xyz)[atom.name.strip().upper()] = atom.xyz
  # }}}

  # {{{ set_sample_points
  def set_sample_points(self, sample_factor = 1./4) :
    sample_spacing = sample_factor*self.resolution
    base = {}
    rn = self.residue.resname.strip().upper()
    if rn in ['DA', 'DG'] :
      base['C2'] = self.xyz.C2
      base['N1'] = self.xyz.N1
      base['C6'] = self.xyz.C6
    elif rn in ['DT', 'DC'] :
      base['N3'] = self.xyz.C2
      base['C4'] = self.xyz.N1
      base['C5'] = self.xyz.C6
    self.sample_points = box_base.get_points_around_base(base,rn,sample_spacing)
  # }}}

  # {{{ get_sample_points_kin_group
  def get_sample_points_kin_group(self,color='red') :
    gn = '%s %s' % (self.chain_id, self.residue.resseq)
    return box_base.get_kin_balls(self.sample_points,
                                  color=color,
                                  group_name=gn)
  # }}}

  # {{{ get_density_sample_points_kin_group
  def get_density_sample_points_kin_group(self,positive_threshold,
                                               negative_threshold) :
    group_name = '%s %s' % (self.chain_id, self.residue.resseq)
    kin = '@group {%s}\n' % group_name
    negative_color = 'orange'
    positive_color = 'sky'
    inb_color = 'white'
    neg_list, pos_list, inbetween = [], [], []
    for k,point in self.density_sample_points.sample_points.items() :
      # print 'point.map_value : %.3f %s' % (point.map_value,  point.map_value > 2)
      id = '%.3f %s %s %s' %(point.map_value,point.plane,point.column,point.row)
      kin_line = '{%s} P %.3f %.3f %.3f\n' % (id,point.xyz[0],
                                          point.xyz[1],point.xyz[2])
      if point.map_value >= positive_threshold : pos_list.append(kin_line)
      elif point.map_value <= negative_threshold : neg_list.append(kin_line)
      else : inbetween.append(kin_line)
    kin += '@balllist {>=%f} color= %s radius = 0.02\n' % (positive_threshold,positive_color)
    for kl in pos_list : kin += kl
    kin += '@balllist {<=%f} color= %s radius = 0.02\n' % (negative_threshold,negative_color)
    for kl in neg_list : kin += kl
    kin += '@balllist {inbetween} color= %s radius = 0.02\n' % inb_color
    for kl in inbetween : kin += kl
    return kin
  # }}}

  # {{{ get_peaks
  def get_peaks(self,positive_threshold, negative_threshold, count_threshold) :
    self.density_sample_points.identify_diff_peaks(positive_threshold, 
                                                   negative_threshold,
                                                   count_threshold)
    # Put peaks in a list that contains peak objects that can give more info 
    # about the peak's location relative to the base atoms
    self.positive_peaks = []
    self.negative_peaks = []
    for peak in self.density_sample_points.positive_peaks :
      self.positive_peaks.append(Peak(peak,self))
    for peak in self.density_sample_points.negative_peaks :
      self.negative_peaks.append(Peak(peak,self))
  # }}}
#}}}

class SimpleBaseClass(object) :

  def __init__(self, chain_id, residue, altloc, resolution) :
    self.chain_id   = chain_id
    self.res_num = residue.resseq
    self.atom_names = [atom.name.strip().upper() for atom in residue.atoms()]
#    self.residue    = residue % can't pickle
    self.alt_loc     = altloc
    self.res_type = residue.resname.strip()
    self.ins_code = residue.icode
    self.resolution = resolution
    self.xyz = group_args()
    self.density_sample_points = DensitySamplePoints()
    self.set_atoms_xyz(residue)
    self.set_sample_points()

  def is_purine(self) :
    if self.res_type.upper() in ['DA','DG'] : return True
    return False

  def is_pyrimidine(self) :
    if self.res_type.upper() in ['DC','DT'] : return True
    return False

  def get_string_id(self,delimiter=',') :
    d = delimiter
    s = '%s%s%s%s%s%s%s%s%s'
    return s % (self.chain_id,d,
                self.res_num,d,
                self.ins_code,d,
                self.alt_loc,d,
                self.res_type)

  def set_atoms_xyz(self, residue, atom_names=None) :
    # atom_names should be a list of string atom names, if None gets all atoms
    # returns a group_args with atom_name->(x,y,z)
    for atom in residue.atoms():
      if atom_names == None :
        vars(self.xyz)[atom.name.strip().upper()] = atom.xyz
      elif atom.name.strip() in atom_names :
        vars(self.xyz)[atom.name.strip().upper()] = atom.xyz

  def set_sample_points(self, sample_factor = 1./4) :
    sample_spacing = sample_factor*self.resolution
    base = {}
#     print self.get_string_id()
    rn = self.res_type.strip().upper()
    if rn in ['DA', 'DG'] :
      base['C2'] = self.xyz.C2
      base['N1'] = self.xyz.N1
      base['C6'] = self.xyz.C6
      # test_points is where we test for positive density when looking for 
      # potential decoys. 
      tp = box_base.get_test_points(base,rn,sample_spacing)
      vars(self.xyz)['test_1'] = tp.point_1
      vars(self.xyz)['test_2'] = tp.point_2
    elif rn in ['DT', 'DC'] :
      base['N3'] = self.xyz.C2
      base['C4'] = self.xyz.N1
      base['C5'] = self.xyz.C6
    self.sample_points = box_base.get_points_around_base(base,rn,sample_spacing)
    

  def get_sample_points_kin_group(self,color='red') :
    gn = '%s %s' % (self.chain_id, self.residue.resseq)
    return box_base.get_kin_balls(self.sample_points,
                                  color=color,
                                  group_name=gn)

  def get_density_sample_points_kin_group(self,positive_threshold,
                                               negative_threshold) :
    group_name = '%s %s' % (self.chain_id, self.res_num)
    kin = '@group {%s}\n' % group_name
    negative_color = 'orange'
    positive_color = 'sky'
    inb_color = 'white'
    neg_list, pos_list, inbetween = [], [], []
    for k,point in self.density_sample_points.sample_points.items() :
      # print 'point.map_value : %.3f %s' % (point.map_value,  point.map_value > 2)
      id = '%.3f %s %s %s' %(point.map_value,point.plane,point.column,point.row)
      kin_line = '{%s} P %.3f %.3f %.3f\n' % (id,point.xyz[0],
                                          point.xyz[1],point.xyz[2])
      if point.map_value >= positive_threshold : pos_list.append(kin_line)
      elif point.map_value <= negative_threshold : neg_list.append(kin_line)
      else : inbetween.append(kin_line)
    kin += '@balllist {>=%f} color= %s radius = 0.02\n' % (positive_threshold,positive_color)
    for kl in pos_list : kin += kl
    kin += '@balllist {<=%f} color= %s radius = 0.02\n' % (negative_threshold,negative_color)
    for kl in neg_list : kin += kl
    kin += '@balllist {inbetween} color= %s radius = 0.02\n' % inb_color
    for kl in inbetween : kin += kl
    return kin

  def get_peaks(self,positive_threshold, negative_threshold, count_threshold) :
    self.density_sample_points.identify_diff_peaks(positive_threshold, 
                                                   negative_threshold,
                                                   count_threshold)
    # Put peaks in a list that contains peak objects that can give more info 
    # about the peak's location relative to the base atoms
    self.positive_peaks = []
    self.negative_peaks = []
    for peak in self.density_sample_points.positive_peaks :
      self.positive_peaks.append(Peak(peak,self))
    for peak in self.density_sample_points.negative_peaks :
      self.negative_peaks.append(Peak(peak,self))

# {{{ DiffDensityAroundBases
class DiffDensityAroundBases(object) :

  # {{{ __init__
  def __init__(self, pdb_id,  positive_threshold = 3,
                              negative_threshold = -3,
                              count_threshold = 4, 
                              log=sys.stderr) :  
    self.log = log
    self.pdb_id = pdb_id
    self.positive_threshold = positive_threshold
    self.negative_threshold = negative_threshold
    self.count_threshold    = count_threshold
    # save model and data in current directory 
    self.mc = GetMapCoeffs(args = [self.pdb_id,'prefix=%s_ddab' % self.pdb_id])
    # the mtz and pdb files are mc.mtz_file mc.pdb_file, respectively
    # get map coefficients
    self.mc.compute_map()
    # print dir(self.mc);sys.exit()
    # the map coefficients file is mc.map_coeff_file_name
    # the resolution range is mc.resolution_range
    # get points around base
    # print self.mc.resolution_range;sys.exit()
    self.resolution = self.mc.resolution_range[1]
    self.set_dna_bases()
    self.get_density_values() # finds mFo-DFc peaks

  def set_thresholds(self, positive = 3, negative = -3, count = 4):
    self.positive_threshold = positive
    self.negative_threshold = negative
    self.count_threshold    = count

  def set_dna_bases(self) :
    self.dna_bases = []
    for chain in self.mc.pdb_in.hierarchy.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            if residue.resname.strip() in ['DA','DT','DC','DG'] :
              b = SimpleBaseClass(chain.id, residue, conformer.altloc,
                self.mc.resolution_range[1] )
              self.dna_bases.append(b)
  # }}}

  # {{{ get_miller_arrays_3dmap
  def get_miller_arrays_3dmap(self) :
    args = [self.mc.map_coeff_file_name,self.mc.pdb_file, 'label=FOFCWT']
    print >> self.log, '='*79
    print >> self.log, 'getting map coeff miller arrays and map'
    print >> self.log, '='*79
    ma,map_3d = map_value_at_point.run(args,sys.stdout)
    return ma,map_3d
  # }}}

  # {{{ get_density_values : gets peaks too
  def get_density_values(self) :
    ma,map_3d = self.get_miller_arrays_3dmap()
    for base in self.dna_bases :
      for k,xyz in base.sample_points.items() :
        point_frac = ma.unit_cell().fractionalize(xyz)
        map_value  = map_3d.eight_point_interpolation(point_frac)
        ga = group_args(plane     = k[0],
                        column    = k[1],
                        row       = k[2],
                        xyz       = xyz,
                        map_value = map_value)
        base.density_sample_points.deposit_sample_point(ga)
      base.get_peaks(self.positive_threshold, 
                     self.negative_threshold,
                     self.count_threshold)
  # }}}

  # {{{ write_all_peak_summary
  # writes a sommary of ALL peaks near a base
  def write_all_peak_summary(self, adjacency_radius=1.5, log=sys.stdout) :
    print >> log, '*'*32 + ' Begin summary ' + '*'*32
    print >> log, 'PARAMETERS'
    print >> log, '  positive_threshold : %s' % self.positive_threshold
    print >> log, '  negative_threshold : %s' % self.negative_threshold
    print >> log, '  count_threshold : %s' % self.count_threshold
    print >> log, '*'*12 + ' BASES ' + '*'*12
    for base in self.dna_bases :
      atoms = base.residue.atoms()
      ppeaks = base.positive_peaks
      npeaks = base.negative_peaks
      if len(ppeaks) == 0 and len(npeaks) == 0 : continue
      print >> log, '  ' + '*'*77
      print >> log, '  %s' % base.get_string_id()
      print >> log, '    Positive peaks : %i' % len(ppeaks)
      if len(ppeaks) > 0 :
        for i,peak in enumerate(ppeaks) :
          print >> log, '      Peak %i # of points : %i' % (i+1, peak.n_sample_points)
        # atom adjacency info
        for atom in atoms :
          a = atom.name.strip().upper()
          pa = [peak for peak in ppeaks if peak.adjacent_to(a, adjacency_radius)]
          print >> log, '        Peaks adjacent to %s : %i ' % (a,len(pa))
      print >> log, '    Negative peaks : %i' % len(npeaks)
      if len(npeaks) > 0 :
        for i,peak in enumerate(npeaks) :
          print >> log, '      Peak %i # of points : %i' % (i+1, peak.n_sample_points)
        # atom adjacency info
        for atom in atoms :
          a = atom.name.strip().upper()
          pa = [peak for peak in npeaks if peak.adjacent_to(a, adjacency_radius)]
          print >> log, '        Peaks adjacent to %s : %i ' % (a,len(pa))
    print >> log, '*'*32 + '  End summary  ' + '*'*32
  # }}}

  # {{{ set_potential_decoys
  def set_potential_decoys(self, 
                           test_radii = 1, 
                           N7_radii=1.9, 
                           C4_radii=1.9, 
                           log=sys.stderr) :
    # This function looks for potential mis-modeled bases by looking for 
    # positive/negative mFo-DFc regions in the following areas : 
    #    - Around test points near N7 and C4 of the purine (positive)
    #    - Around N1, C2, and C8 of the purine (negative)
    # This means that the script is looking for regions of difference density 
    # around 5 different points. 
    # Super strong decoy candidates have a negative peak near C8 and 2
    #   additional peaks.
    # Strong decoy candidates are not super strong that have either 
    #   positive peaks near test 1 AND 2
    #      OR
    #   has at least three of the five ponts with peaks near by.
    # Weak decoys has at least two of the five ponts with peaks near by. 
    # Note that the radii
    # given as arguments are not true radii. It will only look within the sample
    # space around the base. The sample ponts are distributed in planes 
    # parellel to the base and the planes extend at most 0.7 anstroms above and
    # below the base. This means that a radius of 2 will not catch any peaks
    # 0.7 A or more above or below (relatiive to the base plane) the test points.
    print >> log, '*'*36 + ' Decoys ' + '*'*35
    pdbres = 'PDB : %s    Resolution : %s'
    print >> log, pdbres % (self.pdb_id,self.resolution)
    print >> log, 'PARAMETERS'
    print >> log, '  positive_threshold : %s' % self.positive_threshold
    print >> log, '  negative_threshold : %s' % self.negative_threshold
    print >> log, '  count_threshold : %s' % self.count_threshold
#    print >> log, '  N7_radii : %s' % N7_radii
#    print >> log, '  C4_radii : %s' % C4_radii
    print >> log, '  test_radii : %s' % test_radii
    print >> log, '*'*12 + ' BASES ' + '*'*12
    self.weak_decoys = []
    self.strong_decoys = []
    self.super_strong_decoys = []
    # GET DIFFERENCE PEAKS
    for base in self.dna_bases :
      if not base.is_purine() : continue
      ppeaks = base.positive_peaks
      npeaks = base.negative_peaks
      if len(ppeaks) == 0 and len(npeaks) == 0: continue
      # atoms = base.residue.atoms()
      # print >> log, '  ' + '*'*77
      # print >> log, '  %s' % base.get_string_id()
      # print >> log, '    Positive peaks : %i' % len(ppeaks)
#      # get peaks near N7
#       P1p = [p for p in ppeaks if p.adjacent_to('N7', N7_radii)]
#       # get peaks near C4
#       P2p = [p for p in ppeaks if p.adjacent_to('C4', C4_radii)]
      P1p = [p for p in ppeaks if p.adjacent_to('test_1', test_radii)]
      P2p = [p for p in ppeaks if p.adjacent_to('test_2', test_radii)]
      # get negative peaks near C2
      NC2p = [p for p in npeaks if p.adjacent_to('C2', test_radii)]
      # get negative peaks near N1
      NN1p = [p for p in npeaks if p.adjacent_to('N1', test_radii)]
      # get negative peaks near C8
      NC8p = [p for p in npeaks if p.adjacent_to('C8', test_radii)]
      if len(P1p) == 0 and len(P2p) == 0 : continue
      # check to see if there are peaks near any of the test points and atoms.
      # If C8 has a negative peak near by and 2 additional peaks near the
      # remaining 4 points deposit it in super_strong_decoys 
      # If there is a positive peak near test 1 and 2 OR three of the five 
      #   test points have a peak near by then deposit it in strong_decoys. 
      # If only two of the five points have a peak near by then deposit it
      #   into weak_decoys.
      ga = group_args(base     = base,
                    P1peaks  = P1p,
                    P2peaks  = P2p,
                    NC2peaks = NC2p,
                    NN1peaks = NN1p,
                    NC8peaks = NC8p)
      i = 0
      for e in [P1p,P2p,NC2p,NN1p,NC8p] :
        if len(e) > 0 : i += 1
      if len(NC8p) > 0 and i >= 3 :
        print >> log, '%s, super_strong_decoy_candidate' % base.get_string_id()
        self.super_strong_decoys.append(ga)
      elif (len(P1p) > 0 and len(P2p) > 0) or (i > 2) :
        print >> log, '%s, strong_decoy_candidate' % base.get_string_id()
        self.strong_decoys.append(ga)
      elif i > 1 :
        print >> log, '%s, weak_decoy_candidate' % base.get_string_id()
        self.weak_decoys.append(ga)
    print >> log, '  End Decoy summary  '.center(79,'*')
  # }}}

  def write_potential_decoy_summaryOLD(self, log=sys.stderr) :
    print >> log, ' Decoys Candidates '.center(79,'*')
    pdbres = 'PDB : %s    Resolution : %s'
    print >> log, pdbres % (self.pdb_id,self.mc.resolution_range[1])
    print >> log, 'PARAMETERS'
    print >> log, '  positive_threshold : %s' % self.positive_threshold
    print >> log, '  negative_threshold : %s' % self.negative_threshold
    print >> log, '  count_threshold : %s' % self.count_threshold
    print >> log, ' BASES '.center(40)
    print >> log, ' # of strong decoy candidates : %i ' % len(self.strong_decoys)
    for base in self.strong_decoys :
      print >> log, '  %s' % base.get_string_id()
    print >> log, ' # of weak decoy candidates : %i ' % len(self.weak_decoys)
    for base in self.weak_decoys :
      print >> log, '  %s' % base.get_string_id()
    print >> log, ' End Decoy summary '.center(79,'*')

  def write_simpe_summary(self, ga, log) :
    pmsg = '   %s positive peak sample points n : %i'
    nmsg = '   %s negative peak sample points n : %i'
    print >> log, '  %s' % ga.base.get_string_id()
    for P1p in ga.P1peaks :
      print >> log, pmsg % ('test 1', P1p.n_sample_points)
    for P2p in ga.P2peaks :
      print >> log, pmsg % ('test 2', P2p.n_sample_points)
    for NC2p in ga.NC2peaks :
      print >> log, nmsg % ('C2', NC2p.n_sample_points)
    for NN1p in ga.NN1peaks :
      print >> log, nmsg % ('N1', NN1p.n_sample_points)
    for NC8p in ga.NC8peaks :
      print >> log, nmsg % ('C8', NC8p.n_sample_points)

  def write_potential_decoy_summary(self, log=sys.stderr) :
    print >> log, ' Decoys Candidates '.center(79,'*')
    pdbres = 'PDB : %s    Resolution : %s'
    print >> log, pdbres % (self.pdb_id,self.resolution)
    print >> log, 'PARAMETERS'
    print >> log, '  positive_threshold : %s' % self.positive_threshold
    print >> log, '  negative_threshold : %s' % self.negative_threshold
    print >> log, '  count_threshold : %s' % self.count_threshold
    print >> log, ' BASES '.center(40,'*').center(79)
    pmsg = '   %s positive peak sample points n : %i'
    nmsg = '   %s negative peak sample points n : %i'
    ns = '\n # of %s decoy candidates : %i '

    print >> log, ns % ('super strong', len(self.super_strong_decoys))
    for ga in self.super_strong_decoys :
      self.write_simpe_summary(ga,log)

    print >> log, ns % ('strong', len(self.strong_decoys))
    for ga in self.strong_decoys :
      self.write_simpe_summary(ga,log)

#       for C4p in ga.C4peaks :
#         print >> log, '   C4 peak sample points n : %i' % (C4p.n_sample_points)
#       for N7p in ga.N7peaks :
#        print >> log, '   N7 peak sample points n : %i' % (N7p.n_sample_points)

    print >> log, ns % ('weak', len(self.weak_decoys))
    for ga in self.weak_decoys :
      self.write_simpe_summary(ga,log)

    print >> log, ' End Decoy summary '.center(79,'*')

  # {{{ clean up files
  def clean_up_files(self) :
    assert len(self.pdb_id) == 4 and self.pdb_id.isalnum()
    libtbx.easy_run.fully_buffered('rm %s*' % self.pdb_id)
  # }}}

# }}}

# {{{ get_sample_points_kin
def get_sample_points_kin(bases,positive_threshold,negative_threshold) :
  kin = ''
  for base in bases :
    #print base.density_sample_points.sample_points;sys.exit()
    num = base.res_num.strip()
    cid = base.chain_id
    if (num == '8' and cid == 'C'):# or (num == '8' and cid == 'C') :
      print base.get_string_id() + ':' + str(base.xyz.test_1)
      print base.get_string_id() + ':' + str(base.xyz.test_2) 
      kin += base.get_density_sample_points_kin_group(positive_threshold,
                                                      # negative_threshold)
    # kin += base.get_density_sample_points_kin_group(positive_threshold,
                                                      negative_threshold)

  return kin
#}}}

def run(args) :
  
  # {{{ parse args
  if len(args) == 0 : usage()
  pdb_id = args[0]
  assert len(pdb_id) == 4
  threshold = 3
  positive_threshold = threshold
  negative_threshold = -1*threshold
  point_n_threshold  = 4
  write_sample_points_kin = False
  for arg in args[1:] :
    if arg.startswith('positive_threshold') : 
      positive_threshold = float(arg.split('=')[1])
    elif arg.startswith('negative_threshold') : 
      negative_threshold = float(arg.split('=')[1])
    elif arg.startswith('threshold') : 
      negative_threshold = float('-%s' % arg.split('=')[1])
      positive_threshold = float(arg.split('=')[1])
    elif arg.startswith('point_n_threshold') :
      point_n_threshold = int(arg.split('=')[1])
    elif arg.startswith('write_sample_points_kin') :
      if arg.split('=')[1].lower() == 'true' :
        write_sample_points_kin = True
  # }}}

  ddab = DiffDensityAroundBases(pdb_id, positive_threshold,
                                        negative_threshold,
                                        count_threshold = point_n_threshold)
  if write_sample_points_kin :
    fn = 'sample_points_%s.kin' % pdb_id
    fle = open(fn,'w')
    fle.write(get_sample_points_kin(ddab.dna_bases, positive_threshold,
                                              negative_threshold))
    fle.close()
  # write summary of diffence poits to stdout
  ddab.set_potential_decoys()
  ddab.clean_up_files()
  ddab.write_potential_decoy_summary()
    

if __name__ == '__main__' :
  run(sys.argv[1:])

