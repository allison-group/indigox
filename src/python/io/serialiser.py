## \file serialiser.py
from pathlib import Path
import indigox as ix

__all__ = ["SaveITPFile"]

def SaveITPFile(path, mol, pmol=None):
  h_mass = ix.GetPeriodicTable()["H"].GetAtomicMass()
  
  # atom printing
  atm_fmt_str = "{index:>5} {type:>6} {residx:>5} {resname:>8} {name:>7} {chargegroup:>5} {charge:>11.5f} {mass:>9.4f}  {extra}"
  for atom in mol.GetAtoms():
    atm_dat = {"index": atom.GetIndex() + 1,
              "type": atom.GetType().GetName() if atom.HasType() else "%%%",
              "residx": 1,
              "resname": "TEST",
              "name": atom.GetName(),
              "chargegroup": atom.GetTag(),
              "charge": atom.GetPartialCharge(),
              "mass": atom.GetElement().GetAtomicMass() + atom.GetImplicitCount() * h_mass
              }
    if pmol is not None:
      patom = pmol.GetAtom(atom)
      extra_info = ";"
      if not len(patom.GetMappedCharges()):
        extra_info += " UNMAPPED"
      else:
        mu = patom.MeanCharge()
        eta = patom.MedianCharge()
        sigma = patom.StandardDeviationCharge()
        if round(mu, 5) != round(eta, 5) or round(sigma, 5) != 0.0:
          extra_info += "mean: {:.5f}, median: {:.5f}, stdev: {:.5f}".format(mu, eta, sigma)
      if extra_info == ";":
        atm_dat["extra"] = ""
      else:
        atm_dat["extra"] = extra_info
    print(atm_fmt_str.format(**atm_dat))

  # bond printing
  bnd_fmt_str = "{atoma:>5} {atomb:>5} {typecode:>5}    gb_{typeid}   {extra}"
  for bond in mol.GetBonds():
    bnd_dat = {"atoma" : bond.GetSourceAtom().GetIndex() + 1,
               "atomb" : bond.GetTargetAtom().GetIndex() + 1,
               "typecode" : 2,
               "typeid" : bond.GetType().GetID() if bond.HasType() else "UNMAPPED"
              }
    if pmol is not None and bond.HasType():
      pbond = pmol.GetBond(bond)
      other = ["gb_{}".format(tm.GetID()) for tm in pbond.GetMappedTypeCounts()]
      bnd_dat["extra"] = ", ".join(x for x in other if x != "gb_{}".format(bnd_dat["typeid"]))
    else:
      bnd_dat["extra"] = ""
    if bnd_dat["extra"]:
      bnd_dat["extra"] = "; Other terms: " + bnd_dat["extra"]
    print(bnd_fmt_str.format(**bnd_dat))

  # angle printing
  ang_fmt_str = "{atoma:>5} {atomb:>5} {atomc:>5} {typecode:>5}    ga_{typeid} "
  if pmol is not None:
    ang_fmt_str += "   {extra}"
  for angle in mol.GetAngles():
    atoms = angle.GetAtoms()
    ang_dat = {"atoma" : atoms[0].GetIndex() + 1,
               "atomb" : atoms[1].GetIndex() + 1,
               "atomc" : atoms[2].GetIndex() + 1,
               "typecode" : 2,
               "typeid" : angle.GetType().GetID() if angle.HasType() else "UNMAPPED",
               "extra" : ""
               }
    if pmol is not None and angle.HasType():
      pangle = pmol.GetAngle(angle)
      other = ["ga_{}".format(tm.GetID()) for tm in pangle.GetMappedTypeCounts()]
      ang_dat["extra"] = ", ".join(x for x in other if x != "ga_{}".format(ang_dat["typeid"]))
    if ang_dat["extra"]:
      ang_dat["extra"] = "; Other terms: " + ang_dat["extra"]
    print(ang_fmt_str.format(**ang_dat))


