from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import json

class ConformationSearcher:
    def __init__(self):
        pass
    
    def detect_large_rings(self, mol, min_size=10):
        ri = mol.GetRingInfo()
        atom_rings = ri.AtomRings()
        large_rings = [ring for ring in atom_rings if len(ring) >= min_size]
        return len(large_rings) > 0, len(large_rings), max([len(ring) for ring in atom_rings], default=0)
    
    def generate_3d_coordinates(self, mol, num_confs=10, method='etkdg'):
        mol = Chem.AddHs(mol)
        
        has_large_ring, num_large_rings, max_ring_size = self.detect_large_rings(mol)
        
        if method == 'etkdg':
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.pruneRmsThresh = 0.5
            params.maxAttempts = 1000
            params.numThreads = 0
            
            if has_large_ring:
                params.useMacrocycleTorsions = True
                params.useSmallRingTorsions = True
                params.useBasicKnowledge = True
                params.useExpTorsionAnglePrefs = True
                params.maxAttempts = 5000
                params.pruneRmsThresh = 1.0
                
                if max_ring_size >= 15:
                    params.useRandomCoords = True
                    params.maxAttempts = 10000
                    params.pruneRmsThresh = -1.0
            
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
            
            if len(cids) == 0 and has_large_ring:
                params.useRandomCoords = True
                params.maxAttempts = 20000
                params.pruneRmsThresh = -1.0
                params.enforceChirality = False
                cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
            
        elif method == 'distance_geometry':
            params = AllChem.DistanceGeometry()
            params.randomSeed = 42
            params.maxAttempts = 1000
            params.numThreads = 0
            
            if has_large_ring:
                params.maxAttempts = 5000
                params.useRandomCoords = True
            
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return mol, cids
    
    def calculate_energy(self, mol, conf_id):
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
            if ff is not None:
                energy = ff.CalcEnergy()
                return energy
        except:
            pass
        return None
    
    def optimize_conformation(self, mol, conf_id, max_iterations=1000):
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
            if ff is not None:
                ff.Minimize(maxIts=max_iterations)
                energy = ff.CalcEnergy()
                return mol, energy
        except:
            pass
        return mol, None
    
    def conformation_to_dict(self, mol, conf_id, energy=None, optimize=True):
        if optimize:
            mol, energy = self.optimize_conformation(mol, conf_id)
        
        conf = mol.GetConformer(conf_id)
        atoms = []
        bonds = []
        
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            atoms.append({
                'index': i,
                'element': atom.GetSymbol(),
                'x': float(pos.x),
                'y': float(pos.y),
                'z': float(pos.z),
                'atomic_num': atom.GetAtomicNum()
            })
        
        for bond in mol.GetBonds():
            bonds.append({
                'begin': bond.GetBeginAtomIdx(),
                'end': bond.GetEndAtomIdx(),
                'order': int(bond.GetBondType())
            })
        
        return {
            'atoms': atoms,
            'bonds': bonds,
            'energy': energy,
            'conf_id': conf_id,
            'smiles': Chem.MolToSmiles(mol)
        }
    
    def search_conformations(self, smiles, method='etkdg', num_confs=10, callback=None):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        has_large_ring, num_large_rings, max_ring_size = self.detect_large_rings(mol)
        
        mol, cids = self.generate_3d_coordinates(mol, num_confs=num_confs, method=method)
        
        results = []
        valid_confs = [cid for cid in cids if cid >= 0]
        
        if len(valid_confs) == 0:
            ring_info = f" (检测到{num_large_rings}个大环，最大环大小: {max_ring_size}元环)" if has_large_ring else ""
            raise RuntimeError(f"No conformers generated{ring_info}. Try adjusting parameters or using a different method.")
        
        for i, cid in enumerate(valid_confs):
            energy = self.calculate_energy(mol, cid)
            conf_dict = self.conformation_to_dict(mol, cid, energy=energy, optimize=True)
            
            results.append(conf_dict)
            
            if callback:
                callback(conf_dict, i, len(valid_confs))
        
        results.sort(key=lambda x: x['energy'] if x['energy'] is not None else float('inf'))
        
        for i, result in enumerate(results):
            result['rank'] = i + 1
        
        return results

if __name__ == '__main__':
    searcher = ConformationSearcher()
    results = searcher.search_conformations('CCO', num_confs=5)
    for r in results:
        print(f"Rank {r['rank']}: Energy = {r['energy']:.2f} kJ/mol")
