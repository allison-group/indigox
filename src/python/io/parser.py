## \file parser.py
from pathlib import Path
import indigox as ix

__all__ = ["LoadGromosInteractionFunctionParameterFile",
           "LoadPDBFile", "LoadMTBFile", "LoadParameterisedMolecule"]

## \brief Loads a GROMOS IFP file as a forcefield
#  \details See the GROMOS Manual for definition of the IFP file format.
#  \param path the path to the IFP file
#  \return the loaded forcefield
def LoadGromosInteractionFunctionParameterFile(path):
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
        ff.NewAtomType(int(dat[0]), dat[1])
  return ff

## \brief Loads the given PDB file into a Molecule
#  \return the loaded molecule
def LoadPDBFile(path):
  mol = ix.CreateMolecule()
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
        
  return mol

def LoadMTBFile(path, ff = None):
  if ff is None:
    ff = ix.GenerateGROMOS54A7()
  file = list(ix.LoadFile(path.expanduser()))

  blocks = [file[0]]
  for i in range(len(file) - 1):
    if file[i] == "END" and len(file[i+1]):
      blocks.append(file[i+1])
    
  mol = ix.CreateMolecule()
  bondtype = ix.BondType.Harmonic
  angletype = ix.AngleType.Harmonic
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
          atom.SetImplicitCount(atom.GetType().GetImplicitHydrogenCount())
        elif i < start_line[1] + counts[1]:  # Bond data
          bond = mol.NewBond(mol.GetAtomTag(int(dat[0])),
                             mol.GetAtomTag(int(dat[1])))
          bond.SetType(ff.GetBondType(bondtype, int(dat[2])))
        elif i < start_line[2] + counts[2]:  # Angle data
          angle = mol.GetAngle(mol.GetAtomTag(int(dat[0])),
                               mol.GetAtomTag(int(dat[1])),
                               mol.GetAtomTag(int(dat[2])));
          angle.SetType(ff.GetAngleType(angletype, int(dat[2])))
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
          dihedral.AddType(ff.GetDihedralType(propertype, int(dat[4])))

  return mol


def LoadParameterisedMolecule(coord_path, param_path, ff):
  if set(coord_path.suffixes) not in [{".pdb"}, {".aa",".pdb"}, {".ua",".pdb"}]:
    raise FileNotFoundError("Unable to handle file: {}".format(coord_path.name))
  if set(param_path.suffixes) not in [{".mtb"}, {".aa",".mtb"}, {".ua",".mtb"}]:
    raise FileNotFoundError("Unable to handle file: {}".format(param_path.name))
  # assume all atoms the same
  mol_coords = LoadPDBFile(coord_path) # only gives coordinate information
  mol_params = LoadMTBFile(param_path, ff)
  if mol_coords.NumAtoms() != mol_params.NumAtoms():
    raise TypeError("Input files have mismatch atom counts")
  for atom in mol_params.GetAtoms():
    c_atom = mol_coords.GetAtomTag(atom.GetTag())
    atom.SetX(c_atom.GetX())
    atom.SetY(c_atom.GetY())
    atom.SetZ(c_atom.GetZ())

  return mol_params

## \cond
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
