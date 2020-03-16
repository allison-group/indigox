## \file parser.py
from pathlib import Path
from collections import defaultdict
import indigox as ix
import os

__all__ = ["LoadIFPFile", "LoadPDBFile", "LoadMTBFile", "LoadIXDFile",
           "LoadParameterisedMolecule", "LoadITPFile" , "LoadFragmentFile",
           "LoadSDFFile"]

## \brief Loads a GROMOS IFP file as a forcefield
#  \details See the GROMOS Manual for definition of the IFP file format.
#  \param path the path to the IFP file
#  \return the loaded forcefield
def LoadIFPFile(path):
  print("Loading IFP file from path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  pt = ix.GetPeriodicTable()
  data = list(ix.LoadFile(path.expanduser()))
  blocks = [data[0]]
  for i in range(len(data) - 1):
    if data[i] == 'END':
      blocks.append(data[i+1])
  ff = ix.Forcefield(ix.FFFamily.GROMOS, path.stem)
  for b in blocks:
    start = data.index(b) + 1
    end = data[start:].index('END') + start
    if b == 'BONDSTRETCHTYPECODE':
      ff.ReserveBondTypes(ix.BondType.Harmonic,
                          int(data[start].split()[0]))
      ff.ReserveBondTypes(ix.BondType.Quartic,
                          int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, kq, kh, b0 = map(float, line.split())
        ff.LinkBondTypes(ff.NewBondType(ix.BondType.Harmonic, int(idx), kh, b0),
                         ff.NewBondType(ix.BondType.Quartic, int(idx), kq, b0))
    elif b == 'BONDANGLEBENDTYPECODE':
      ff.ReserveAngleTypes(ix.AngleType.Harmonic,
                           int(data[start].split()[0]))
      ff.ReserveAngleTypes(ix.AngleType.CosineHarmonic,
                           int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, kch, kh, t0 = map(float, line.split())
        ff.LinkAngleTypes(ff.NewAngleType(ix.AngleType.Harmonic, int(idx), kh, t0),
                          ff.NewAngleType(ix.AngleType.CosineHarmonic, int(idx), kch, t0))
    elif b == 'IMPDIHEDRALTYPECODE':
      ff.ReserveDihedralTypes(ix.DihedralType.Improper,
                              int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, k, epsilon = map(float, line.split())
        ff.NewDihedralType(ix.DihedralType.Improper, int(idx), k, epsilon)
    elif b == 'TORSDIHEDRALTYPECODE':
      ff.ReserveDihedralTypes(ix.DihedralType.Proper,
                              int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, k, phase, m = map(float, line.split())
        ff.NewDihedralType(ix.DihedralType.Proper, int(idx), k, phase, int(m))
    elif b == 'SINGLEATOMLJPAIR':
      num_lines = int((end - start - 1) / int(data[start]))
      atm_count = int(data[start])
      ff.ReserveAtomTypes(atm_count)
      for i in range(atm_count):
        lines = []
        for j in range(num_lines):
          lines.append(data[start+1+j+i*num_lines])
        dat = "   ".join(x for x in lines).split()
        atm_dat = {'int_code':int(dat[0]),
                   'name':dat[1],
                   'c6':float(dat[2])**2,
                   'c12':[float(dat[3])**2,float(dat[4])**2,float(dat[5])**2,],
                   'c6_1_4':float(dat[6])**2,
                   'c12_1_4':float(dat[7])**2,
                   'interactions':dict(),
                  }
        ff.NewAtomType(int(dat[0]), dat[1], pt[dat[1][0]])
  return ff

## \brief Loads the given PDB file into a Molecule
#  \return the loaded molecule
def LoadPDBFile(path, details = None):
  print("Loading molecule from PDB file at path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  mol = ix.Molecule(path.stem)
  got_atoms = False
  for line in ix.LoadFile(path.expanduser()):
    record_type = line[0:6]
    if record_type in ["HETATM", "ATOM  "] and not got_atoms:
      atom = mol.NewAtom()
      try:
        atom.SetTag(int(line[6:11]))
      except TypeError:
        pass
      atom.SetName(line[12:16].strip())
      atom.SetElement(line[76:78].strip())
      try:
        atom.SetPosition(float(line[30:38])/10,
                         float(line[38:46])/10,
                         float(line[46:54])/10)
      except TypeError:
        atom.SetPosition(0.0,0.0,0.0)
    elif record_type == "CONECT":
      try:
        a = int(line[6:11])
      except ValueError:
        a = None
      try:
        b1 = int(line[11:16])
      except ValueError:
        b1 = None
      try:
        b2 = int(line[16:21])
      except ValueError:
        b2 = None
      try:
        b3 = int(line[21:26])
      except ValueError:
        b3 = None
      try:
        b4 = int(line[26:31])
      except ValueError:
        b4 = None

      if a is None:
        continue
      else:
        try:
          a = mol.GetAtomTag(a)
        except IndexError:
          continue
      for b in [b1,b2,b3,b4]:
        if b is None:
          continue
        else:
          try:
            b = mol.GetAtomTag(b)
          except IndexError:
            continue
        if a.GetIndex() == b.GetIndex():
          continue
        try:
          mol.NewBond(a, b)
        except IndexError:
          continue
    elif record_type == "ENDMDL":
      got_atoms = True

  if details is not None:
    LoadIXDFile(details, mol)
  return mol

def LoadIXDFile(path, mol):
  print("Loading formal charges and bond orders from IXD file at path %s." % str(os.path.join(os.getcwd(), path)))
  set_implicits = []
  path = Path(path)
  for line in ix.LoadFile(path.expanduser()):
    line = line.split()
    l = list(map(int, line[1:]))
    if line[0] == "ATOM":
      atom = mol.GetAtomTag(l[0])
      if not atom:
        raise IndexError("Unknown atom tag: {}".format(l[0]))
      atom.SetFormalCharge(l[1])
      atom.SetImplicitCount(l[2])
      set_implicits.append(l[0])
    elif line[0] == "BOND":
      bond = mol.GetBond(mol.GetAtomTag(l[0]), mol.GetAtomTag(l[1]))
      if not bond:
        raise IndexError("Unknown bond: {}-{}".format(l[0], l[1]))
      if l[2] == 1:
        bond.SetOrder(ix.BondOrder.Single)
      elif l[2] == 2:
        bond.SetOrder(ix.BondOrder.Double)
      elif l[2] == 3:
        bond.SetOrder(ix.BondOrder.Triple)
      elif l[2] == 4:
        bond.SetOrder(ix.BondOrder.Quadruple)
      elif l[2] == 5:
        bond.SetOrder(ix.BondOrder.Aromatic)
      elif l[2] == 6:
        bond.SetOrder(ix.BondOrder.OneAndAHalf)
      elif l[2] == 7:
        bond.SetOrder(ix.BondOrder.TwoAndAHalf)
    elif line[0] == "MOLECULE":
      mol.SetMolecularCharge(l[0])

  set_implicits = set(set_implicits)
  for atom in mol.GetAtoms():
    if atom.GetElement() != "C":
      continue
    if atom.GetTag() in set_implicits:
      continue

    bonded_count = 0
    aromatic_bonds = False
    for bond in atom.GetBonds():
      if bond.GetOrder() == ix.BondOrder.Single:
        bonded_count += 1
      elif bond.GetOrder() == ix.BondOrder.Double:
        bonded_count += 2
      elif bond.GetOrder() == ix.BondOrder.Triple:
        bonded_count += 3
      elif bond.GetOrder() == ix.BondOrder.Quadruple:
        bonded_count += 4
      elif bond.GetOrder() == ix.BondOrder.Aromatic:
        if not aromatic_bonds:
          bonded_count += 1
          aromatic_bonds = True
        bonded_count += 1
      elif bond.GetOrder() == ix.BondOrder.OneAndAHalf:
        bonded_count += 1.5
      elif bond.GetOrder() == ix.BondOrder.TwoAndAHalf:
        bonded_count += 2.5

    if (((bonded_count - int(bonded_count)) != 0) or bonded_count > 4):
      raise ValueError("Atom tag {} has weird valence state.".format(atom.GetTag()))
    atom.SetImplicitCount(4 - bonded_count)

  tot_charge = sum(atm.GetFormalCharge() for atm in mol.GetAtoms())
  if tot_charge != mol.GetMolecularCharge():
    raise ValueError("Sum of atom formal charges does not match molecular charge")

def LoadMTBFile(path, ff, details=None):
  print("Loading MTB file at path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  file = list(ix.LoadFile(path.expanduser()))

  blocks = [file[0]]
  for i in range(len(file) - 1):
    if file[i] == "END" and len(file[i+1]):
      blocks.append(file[i+1])
    
  mol = ix.Molecule(path.stem)
  bondtype = ix.BondType.Quartic
  angletype = ix.AngleType.CosineHarmonic
  impropertype = ix.DihedralType.Improper
  propertype = ix.DihedralType.Proper
  for b in blocks:
    start = file.index(b) + 1
    
    if b == "MTBUILDBLSOLUTE":
      counts = [0, 0, 0, 0, 0, 0]
      start_line = [0, 0, 0, 0, 0, 0]
      mol.SetName(file[start])
      for i in range(len(counts)):
        pre = start + 1 + len(counts[:i]) + sum(counts[:i])
        counts[i] = int(file[pre].split()[0])
        start_line[i] = pre + 1

      for i in range(len(file)):
        dat = file[i].split()
        if i in (x - 1 for x in start_line) or i <= start:
          continue
        elif i < start_line[0] + counts[0]:  # Atom data
          atom = mol.NewAtom()
          atom.SetTag(int(dat[0]))
          atom.SetName(dat[1])
          atom.SetType(ff.GetAtomType(int(dat[2])))
          atom.SetPartialCharge(float(dat[4]))
          atom.SetElement(atom.GetType().GetElement())
        elif i < start_line[1] + counts[1]:  # Bond data
          bond = mol.NewBond(mol.GetAtomTag(int(dat[0])),
                             mol.GetAtomTag(int(dat[1])))
          bond.SetType(ff.GetBondType(bondtype, int(dat[2])))
        elif i < start_line[2] + counts[2]:  # Angle data
          angle = mol.GetAngle(mol.GetAtomTag(int(dat[0])),
                               mol.GetAtomTag(int(dat[1])),
                               mol.GetAtomTag(int(dat[2])));
          if not angle:
            raise IndexError("No angle found: {}-{}-{}".format(*dat[:3]))
          angle.SetType(ff.GetAngleType(angletype, int(dat[3])))
        elif i < start_line[3] + counts[3]:  # Improper dihedrals
          a = mol.GetAtomTag(int(dat[0]))
          b = mol.GetAtomTag(int(dat[1]))
          c = mol.GetAtomTag(int(dat[2]))
          d = mol.GetAtomTag(int(dat[3]))
          if mol.HasDihedral(a, b, c, d):
            dihedral = mol.GetDihedral(a, b, c, d)
          else:
            dihedral = mol.NewDihedral(a, b, c, d)
          dihedral.AddType(ff.GetDihedralType(impropertype, int(dat[4])))
        elif i < start_line[4] + counts[4]:  # Proper dihedrals
          a = mol.GetAtomTag(int(dat[0]))
          b = mol.GetAtomTag(int(dat[1]))
          c = mol.GetAtomTag(int(dat[2]))
          d = mol.GetAtomTag(int(dat[3]))
          dihedral = mol.GetDihedral(a, b, c, d)
          if not dihedral:
            raise IndexError("No dihedral found: {}-{}-{}-{}".format(*dat[:4]))
          dihedral.AddType(ff.GetDihedralType(propertype, int(dat[4])))
  
  if details is not None:
    LoadIXDFile(details, mol)
  
  return mol


def LoadParameterisedMolecule(coord_path, param_path, ff, details=None):
  coord_path = Path(coord_path)
  param_path = Path(param_path)
  if set(coord_path.suffixes) not in [{".pdb"}, {".aa",".pdb"}, {".ua",".pdb"}]:
    raise FileNotFoundError("Unable to handle file: {}".format(coord_path.name))
  if set(param_path.suffixes) not in [{".mtb"}, {".aa",".mtb"}, {".ua",".mtb"},
                                      {".itp"}, {".aa",".itp"}, {".ua",".itp"},
                                      {".top"}]:
    raise FileNotFoundError("Unable to handle file: {}".format(param_path.name))
  # assume all atoms the same
  mol_coords = LoadPDBFile(coord_path) # only gives coordinate information
  if str(param_path).endswith(".mtb"):
    mol_params = LoadMTBFile(param_path, ff, details)
  else:
    mol_params = LoadITPFile(param_path, ff, details)
  if mol_coords.NumAtoms() != mol_params.NumAtoms():
    raise TypeError("Input files have mismatch atom counts")
  for atom in mol_params.GetAtoms():
    c_atom = mol_coords.GetAtomTag(atom.GetTag())
    atom.SetPosition(c_atom.GetX(), c_atom.GetY(), c_atom.GetZ())
    atom.SetElement(c_atom.GetElement())
  return mol_params

def LoadITPFile(path, ff, details=None):
  print("Loading ITP file at path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  mol = ix.Molecule(path.stem)
  
  bondtype = ix.BondType.Harmonic
  angletype = ix.AngleType.Harmonic
  impropertype = ix.DihedralType.Improper
  propertype = ix.DihedralType.Proper
  
  current_block = None
  for line in ix.LoadFile(path.expanduser(), comment=";"):
    if "[" in line:
      data = line.split()
      current_block = data[1]
    elif "#" in line:
      current_block = None
    elif current_block == "atoms":
      data = line.split()
      atom = mol.NewAtom()
      atom.SetTag(int(data[0]))
      atom.SetName(data[4])
      atom.SetType(ff.GetAtomType(data[1]))
      atom.SetPartialCharge(float(data[6]))
      atom.SetElement(atom.GetType().GetElement())
    elif current_block == "bonds":
      data = line.split()
      bond = mol.NewBond(mol.GetAtomTag(int(data[0])),
                         mol.GetAtomTag(int(data[1])))
      if data[-1].startswith("gb"):
        bond.SetType(ff.GetBondType(bondtype, int(data[-1].split("_")[1])))
    elif current_block == "angles":
      data = line.split()
      a1 = mol.GetAtomTag(int(data[0]))
      a2 = mol.GetAtomTag(int(data[1]))
      a3 = mol.GetAtomTag(int(data[2]))
      angle = mol.GetAngle(mol.GetAtomTag(int(data[0])),
                           mol.GetAtomTag(int(data[1])),
                           mol.GetAtomTag(int(data[2])))
      if not angle:
        raise IndexError("No angle found: {}-{}-{}".format(*data[:3]))
      if data[-1].startswith("ga"):
        ty = ff.GetAngleType(angletype, int(data[-1].split("_")[1]))
        angle.SetType(ty)
    elif current_block == "dihedrals":
      data = line.split()
      a = mol.GetAtomTag(int(data[0]))
      b = mol.GetAtomTag(int(data[1]))
      c = mol.GetAtomTag(int(data[2]))
      d = mol.GetAtomTag(int(data[3]))
      if mol.HasDihedral(a,b,c,d):
        dihedral = mol.GetDihedral(a,b,c,d)
      else:
        dihedral = mol.NewDihedral(a,b,c,d)
      if data[-1].startswith("gi"):
        dihedral.AddType(ff.GetDihedralType(impropertype, int(data[-1].split("_")[1])))
      elif data[-1].startswith("gd"):
        dihedral.AddType(ff.GetDihedralType(propertype, int(data[-1].split("_")[1])))

  if details is not None:
    LoadIXDFile(details, mol)
  return mol

def LoadFragmentFile(path, ff):
  print("Loading fragment file at path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  file_data = list(ix.LoadFile(path.expanduser()))
  if file_data.count("MOLECULE") != 1:
    raise InputError("Expect only one molecule per file")
  if file_data[0] != "MOLECULE":
    raise InputError("Expected MOLECULE block")
  if file_data.count("FRAGMENT") != file_data.count("OVERLAP"):
    raise InputError("Should have same number of fragment and overlaps")

  if "END" in file_data[1:4]:
    raise InputError("Missing input files")
  coord_path = path.parent / file_data[1].strip()
  param_path = path.parent / file_data[2].strip()
  detail_path = path.parent / file_data[3].strip()
  
  mol = LoadParameterisedMolecule(coord_path, param_path, ff, detail_path)

  fragments = []
  current_block = None
  frag = []
  overlap = []
  for line in file_data:
    if line == "END":
      if current_block == "OVERLAP":
        frag = ix.VecAtom(frag)
        overlap = ix.VecAtom(overlap)
        fragments.append(ix.Fragment(mol, frag, overlap))
        frag = []
        overlap = []
      current_block = None
    elif line in ["MOLECULE", "FRAGMENT", "OVERLAP"]:
      current_block = line
    elif current_block == "FRAGMENT":
      frag.extend(mol.GetAtomTag(x) for x in map(int, line.split()))
    elif current_block == "OVERLAP":
      overlap.extend(mol.GetAtomTag(x) for x in map(int, line.split()))

  return mol, fragments


## \cond

def LoadSDFFile(path):
  print("Loading SDF file at path %s." % str(os.path.join(os.getcwd(), path)))
  path = Path(path)
  bond_orders = {1: ix.BondOrder.Single,
                 2: ix.BondOrder.Double,
                 3: ix.BondOrder.Triple}
  molecules = []
  pt = ix.GetPeriodicTable()
  data = list(ix.LoadFile(path.expanduser(), comment=None, blanks=True))
  i = 0
  new_molecule = True
  num_atoms = 0
  num_bonds = 0
  while i < len(data):
    if data[i] == "$$$$": 
      i += 1
      new_molecule = True
    elif new_molecule:
      molecules.append(ix.Molecule(data[i]))
      counts = data[i+3].split()
      num_atoms = int(counts[0])
      num_bonds = int(counts[1])
      i += 4
      new_molecule = False
    elif num_atoms:
      atom_line = data[i].split()
      x, y, z = map(float, atom_line[0:3])
      element = pt[atom_line[3]]
      fc = int(atom_line[4])
      atm = molecules[-1].NewAtom(element, x, y, z)
      atm.SetName("{}{}".format(element.GetSymbol(), molecules[-1].NumAtoms()))
      atm.SetFormalCharge(fc)
      atm.SetTag(molecules[-1].NumAtoms())
      num_atoms -= 1
      i += 1
    elif num_bonds:
      bond_line = data[i].split()
      a, b, order, stero = map(int, bond_line)
      atma = molecules[-1].GetAtomTag(a)
      atmb = molecules[-1].GetAtomTag(b)
      molecules[-1].NewBond(atma, atmb).SetOrder(bond_orders[order])
      i+= 1
      num_bonds -= 1
    else:
      i += 1
      
  return molecules

if __name__ == "__main__":
  p = Path('~/tmp/itchy-wookie-data/Input/ForceFields/54A7_original.ifp')
  ff = LoadGromosInteractionFunctionParameterFile(p.expanduser())
  ff2 = ix.GenerateGROMOS54A7()
  print(ff.GetFamily(), ff.GetName())
  print("Atom types:    ", ff.NumAtomTypes())
  print("Bond types:    ", ff.NumBondTypes())
  print("Angle types:   ", ff.NumAngleTypes())
  print("Dihedral types:", ff.NumDihedralTypes())
  print(ff2.GetFamily(), ff2.GetName())
  print("Atom types:    ", ff2.NumAtomTypes())
  print("Bond types:    ", ff2.NumBondTypes())
  print("Angle types:   ", ff2.NumAngleTypes())
  print("Dihedral types:", ff2.NumDihedralTypes())
## \endcond
