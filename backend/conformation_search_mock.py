import math
import random
import time

class MockConformationSearcher:
    def __init__(self):
        random.seed(42)
        
        self.element_data = {
            'C': {'radius': 0.77, 'color': '#333333', 'bonds': 4},
            'H': {'radius': 0.37, 'color': '#FFFFFF', 'bonds': 1},
            'O': {'radius': 0.73, 'color': '#FF0000', 'bonds': 2},
            'N': {'radius': 0.75, 'color': '#0000FF', 'bonds': 3},
            'S': {'radius': 1.02, 'color': '#FFFF00', 'bonds': 2},
            'Cl': {'radius': 0.99, 'color': '#00FF00', 'bonds': 1},
            'F': {'radius': 0.71, 'color': '#00FF00', 'bonds': 1},
        }
    
    def detect_large_rings(self, mol_atoms, mol_bonds, min_size=10):
        from collections import defaultdict, deque
        
        adj = defaultdict(list)
        for bond in mol_bonds:
            adj[bond['begin']].append(bond['end'])
            adj[bond['end']].append(bond['begin'])
        
        max_ring_size = 0
        num_large_rings = 0
        
        def find_rings(start):
            rings = []
            stack = [(start, -1, [start])]
            while stack:
                node, parent, path = stack.pop()
                for neighbor in adj.get(node, []):
                    if neighbor == parent:
                        continue
                    if neighbor in path:
                        ring_start = path.index(neighbor)
                        ring = path[ring_start:]
                        if len(ring) >= 3:
                            rings.append(ring)
                    else:
                        stack.append((neighbor, node, path + [neighbor]))
            return rings
        
        seen_rings = set()
        for atom in mol_atoms:
            rings = find_rings(atom['index'])
            for ring in rings:
                ring_key = tuple(sorted(ring))
                if ring_key not in seen_rings:
                    seen_rings.add(ring_key)
                    ring_size = len(ring)
                    if ring_size > max_ring_size:
                        max_ring_size = ring_size
                    if ring_size >= min_size:
                        num_large_rings += 1
        
        return num_large_rings > 0, num_large_rings, max_ring_size
    
    def smiles_to_atoms(self, smiles):
        atoms = []
        bonds = []
        atom_index = 0
        
        i = 0
        prev_atoms = []
        ring_atoms = {}
        
        while i < len(smiles):
            char = smiles[i]
            
            if char.isupper():
                element = char
                if i + 1 < len(smiles) and smiles[i + 1].islower():
                    element += smiles[i + 1]
                    i += 1
                
                if element in self.element_data:
                    atoms.append({
                        'index': atom_index,
                        'element': element,
                        'atomic_num': self._get_atomic_num(element)
                    })
                    prev_atoms.append(atom_index)
                    atom_index += 1
                i += 1
            elif char == '(':
                branch_start = len(prev_atoms)
                i += 1
            elif char == ')':
                if len(prev_atoms) > 1:
                    prev_atoms.pop()
                i += 1
            elif char.isdigit():
                ring_num = char
                if ring_num in ring_atoms:
                    bonds.append({
                        'begin': ring_atoms[ring_num],
                        'end': atom_index - 1,
                        'order': 1
                    })
                    del ring_atoms[ring_num]
                else:
                    ring_atoms[ring_num] = atom_index - 1
                i += 1
            elif char == '=':
                i += 1
                continue
            else:
                i += 1
        
        for j in range(len(atoms) - 1):
            bonds.append({
                'begin': j,
                'end': j + 1,
                'order': 1
            })
        
        return atoms, bonds
    
    def _get_atomic_num(self, element):
        nums = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16, 'Cl': 17}
        return nums.get(element, 6)
    
    def generate_3d_coordinates(self, atoms, bonds, num_confs=10, has_large_ring=False, max_ring_size=0):
        conformations = []
        
        for conf_idx in range(num_confs):
            conf_atoms = []
            angle_offset = conf_idx * (2 * math.pi / num_confs) * 0.3
            
            if has_large_ring:
                ring_radius = max_ring_size * 0.3
                for idx, atom in enumerate(atoms):
                    angle = idx * (2 * math.pi / max(len(atoms), 1)) + angle_offset
                    height_offset = math.sin(idx * 0.5 + conf_idx * 0.2) * 0.5
                    
                    x = ring_radius * math.cos(angle) + random.gauss(0, 0.15)
                    y = ring_radius * math.sin(angle) + random.gauss(0, 0.15)
                    z = height_offset + (idx - len(atoms) / 2) * 0.3 + random.gauss(0, 0.1)
                    
                    conf_atoms.append({
                        'index': atom['index'],
                        'element': atom['element'],
                        'atomic_num': atom['atomic_num'],
                        'x': x,
                        'y': y,
                        'z': z
                    })
            else:
                for idx, atom in enumerate(atoms):
                    angle = idx * (2 * math.pi / max(len(atoms), 1)) + angle_offset
                    radius = 1.5 + (idx % 3) * 0.5
                    
                    x = radius * math.cos(angle) + random.gauss(0, 0.1)
                    y = radius * math.sin(angle) + random.gauss(0, 0.1)
                    z = (idx - len(atoms) / 2) * 0.5 + random.gauss(0, 0.1)
                    
                    conf_atoms.append({
                        'index': atom['index'],
                        'element': atom['element'],
                        'atomic_num': atom['atomic_num'],
                        'x': x,
                        'y': y,
                        'z': z
                    })
            
            conformations.append(conf_atoms)
        
        return conformations
    
    def calculate_energy(self, atoms):
        energy = 0.0
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dx = atoms[i]['x'] - atoms[j]['x']
                dy = atoms[i]['y'] - atoms[j]['y']
                dz = atoms[i]['z'] - atoms[j]['z']
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                if dist < 1.0:
                    energy += (1.0 / (dist + 0.1)) ** 12
                else:
                    energy += 0.1 / (dist ** 6)
        return energy * 100
    
    def search_conformations(self, smiles, method='etkdg', num_confs=10, callback=None):
        atoms, bonds = self.smiles_to_atoms(smiles)
        
        if not atoms:
            atoms = [
                {'index': 0, 'element': 'C', 'atomic_num': 6},
                {'index': 1, 'element': 'C', 'atomic_num': 6},
                {'index': 2, 'element': 'O', 'atomic_num': 8},
            ]
            bonds = [
                {'begin': 0, 'end': 1, 'order': 1},
                {'begin': 1, 'end': 2, 'order': 1},
            ]
        
        has_large_ring, num_large_rings, max_ring_size = self.detect_large_rings(atoms, bonds)
        
        coords_list = self.generate_3d_coordinates(
            atoms, bonds, num_confs, 
            has_large_ring=has_large_ring, 
            max_ring_size=max_ring_size
        )
        
        results = []
        
        if len(coords_list) == 0:
            ring_info = f" (检测到{num_large_rings}个大环，最大环大小: {max_ring_size}元环)" if has_large_ring else ""
            raise RuntimeError(f"No conformers generated{ring_info}. Try adjusting parameters or using a different method.")
        
        for i, coords in enumerate(coords_list):
            sleep_time = 0.5 if has_large_ring else 0.3
            time.sleep(sleep_time)
            
            energy = self.calculate_energy(coords)
            
            conf_dict = {
                'atoms': coords,
                'bonds': bonds,
                'energy': energy,
                'conf_id': i,
                'smiles': smiles,
                'has_large_ring': has_large_ring,
                'max_ring_size': max_ring_size
            }
            
            results.append(conf_dict)
            
            if callback:
                callback(conf_dict, i, len(coords_list))
        
        results.sort(key=lambda x: x['energy'] if x['energy'] is not None else float('inf'))
        
        for i, result in enumerate(results):
            result['rank'] = i + 1
        
        return results


if __name__ == '__main__':
    searcher = MockConformationSearcher()
    results = searcher.search_conformations('CCO', num_confs=5)
    for r in results:
        print(f"Rank {r['rank']}: Energy = {r['energy']:.2f} kJ/mol")
