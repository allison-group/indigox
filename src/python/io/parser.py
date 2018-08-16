from pathlib import Path
import indigox as ix

def LoadGromosInteractionFunctionParameterFile(path):
  data = list(ix.LoadFile(path.expanduser()))
  blocks = [data[0]]
  for i in range(len(data) - 1):
    if data[i] == 'END':
      blocks.append(data[i+1])
      
  ff = ix.Forcefield(ix.ForcefieldFamily.GROMOS, path.stem)
  for b in blocks:
    start = data.index(b) + 1
    end = data[start:].index('END') + start
    if b == 'BONDSTRETCHTYPECODE':
      ff.ReserveBondTypes(ix.BondFunctionType.Harmonic, 
                          int(data[start].split()[0]))
      ff.ReserveBondTypes(ix.BondFunctionType.Quartic, 
                          int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, kq, kh, b0 = map(float, line.split())
        ff.NewHarmonicBondType(int(idx), kh, b0)
        ff.NewQuarticBondType(int(idx), kq, b0)
    elif b == 'BONDANGLEBENDTYPECODE':
      ff.ReserveAngleTypes(ix.AngleFunctionType.Harmonic,
                           int(data[start].split()[0]))
      ff.ReserveAngleTypes(ix.AngleFunctionType.CosineHarmonic,
                           int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, kch, kh, t0 = map(float, line.split())
        ff.NewHarmonicAngleType(int(idx), kh, t0)
        ff.NewCosineHarmonicAngleType(int(idx), kch, t0)
    elif b == 'IMPDIHEDRALTYPECODE':
      ff.ReserveDihedralTypes(ix.DihedralFunctionType.Improper,
                              int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, k, epsilon = map(float, line.split())
        ff.NewImproperDihedralType(int(idx), k, epsilon)
    elif b == 'TORSDIHEDRALTYPECODE':
      ff.ReserveDihedralTypes(ix.DihedralFunctionType.Proper,
                              int(data[start].split()[0]))
      for line in data[start + 1:end]:
        idx, k, phase, m = map(float, line.split())
        ff.NewProperDihedralType(int(idx), k, phase, int(m))
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

if __name__ == "__main__":
  p = Path('~/tmp/itchy-wookie-data/Input/ForceFields/54A7_original.ifp')
  ff = LoadGromosInteractionFunctionParameterFile(p.expanduser())
  print(ff.GetFamily(), ff.GetName())
  print("Atom types:    ", ff.NumAtomTypes())
  print("Bond types:    ", ff.NumBondTypes())
  print("Angle types:   ", ff.NumAngleTypes())
  print("Dihedral types:", ff.NumDihedralTypes())

