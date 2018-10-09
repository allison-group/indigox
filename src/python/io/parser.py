## \file parser.py
from pathlib import Path
import indigox as ix

__all__ = ["LoadGromosInteractionFunctionParameterFile"]

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
