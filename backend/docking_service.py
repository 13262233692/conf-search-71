import math
import random
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class BindingSite:
    center: List[float]
    radius: float
    name: str = "default"
    residues: List[str] = field(default_factory=list)

@dataclass
class ReceptorAtom:
    index: int
    element: str
    x: float
    y: float
    z: float
    atomic_num: int
    residue_name: str = ""
    residue_id: int = 0
    charge: float = 0.0

@dataclass
class DockingPose:
    atoms: List[Dict]
    score: float
    rotation: List[float]
    translation: List[float]
    rank: int = 0
    clash_count: int = 0
    hbond_count: int = 0

VAN_DER_WAALS_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'S': 1.80, 'P': 1.80, 'Cl': 1.75,
    'Br': 1.85, 'I': 1.98, 'Na': 2.27, 'K': 2.75,
    'Ca': 2.31, 'Mg': 1.73, 'Fe': 1.94, 'Zn': 1.39,
}

ATOMIC_CHARGES = {
    'H': 0.1, 'C': -0.1, 'N': -0.4, 'O': -0.5,
    'S': -0.2, 'F': -0.3, 'Cl': -0.2, 'Br': -0.2,
    'P': 0.3, 'Na': 0.6, 'K': 0.6, 'Ca': 1.0,
    'Mg': 1.0, 'Fe': 0.5, 'Zn': 0.5,
}

PRESET_RECEPTORS = {
    "hydrophobic_pocket": {
        "name": "疏水性结合口袋",
        "description": "典型的疏水性蛋白质结合位点",
        "atoms": [],
        "binding_site": {"center": [0.0, 0.0, 0.0], "radius": 8.0},
        "generate": True,
        "pocket_type": "hydrophobic"
    },
    "polar_pocket": {
        "name": "极性结合口袋",
        "description": "含极性残基的结合位点",
        "atoms": [],
        "binding_site": {"center": [0.0, 0.0, 0.0], "radius": 7.0},
        "generate": True,
        "pocket_type": "polar"
    },
    "metal_binding_site": {
        "name": "金属结合位点",
        "description": "含有金属离子的催化位点",
        "atoms": [],
        "binding_site": {"center": [0.0, 0.0, 0.0], "radius": 6.0},
        "generate": True,
        "pocket_type": "metal"
    }
}


def generate_pocket_atoms(pocket_type="hydrophobic", center=(0, 0, 0), num_residues=8):
    atoms = []
    random.seed(42)
    
    residue_templates = {
        "hydrophobic": [
            ("PHE", [('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6)]),
            ("LEU", [('C', 6), ('C', 6), ('C', 6)]),
            ("VAL", [('C', 6), ('C', 6), ('C', 6)]),
            ("ILE", [('C', 6), ('C', 6), ('C', 6)]),
            ("ALA", [('C', 6)]),
            ("TRP", [('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6)]),
        ],
        "polar": [
            ("SER", [('O', 8), ('C', 6)]),
            ("THR", [('O', 8), ('C', 6)]),
            ("TYR", [('O', 8), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6), ('C', 6)]),
            ("ASN", [('O', 8), ('N', 7), ('C', 6)]),
            ("GLN", [('O', 8), ('N', 7), ('C', 6), ('C', 6)]),
            ("ARG", [('N', 7), ('N', 7), ('C', 6), ('C', 6), ('C', 6)]),
        ],
        "metal": [
            ("HIS", [('N', 7), ('C', 6), ('C', 6), ('C', 6), ('C', 6)]),
            ("CYS", [('S', 16), ('C', 6)]),
            ("ASP", [('O', 8), ('O', 8), ('C', 6), ('C', 6)]),
            ("GLU", [('O', 8), ('O', 8), ('C', 6), ('C', 6), ('C', 6)]),
            ("MG", [('Mg', 12)]),
            ("ZN", [('Zn', 30)]),
            ("FE", [('Fe', 26)]),
        ]
    }
    
    templates = residue_templates.get(pocket_type, residue_templates["hydrophobic"])
    idx = 0
    
    for i in range(num_residues):
        res_template = templates[i % len(templates)]
        res_name = res_template[0]
        
        angle = 2 * math.pi * i / num_residues
        radius = 6.0 + random.gauss(0, 0.3)
        res_center_x = center[0] + radius * math.cos(angle)
        res_center_y = center[1] + radius * math.sin(angle)
        res_center_z = center[2] + random.gauss(0, 0.5)
        
        for atom_element, atomic_num in res_template[1]:
            dx = random.gauss(0, 0.8)
            dy = random.gauss(0, 0.8)
            dz = random.gauss(0, 0.8)
            
            charge = ATOMIC_CHARGES.get(atom_element, 0.0)
            if atom_element in ('Mg', 'Zn', 'Fe'):
                charge = float(atomic_num) * 0.04
            
            atoms.append({
                'index': idx,
                'element': atom_element,
                'x': res_center_x + dx,
                'y': res_center_y + dy,
                'z': res_center_z + dz,
                'atomic_num': atomic_num,
                'residue_name': res_name,
                'residue_id': i + 1,
                'charge': charge
            })
            idx += 1
    
    return atoms


def rotation_matrix(axis, angle):
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c = math.cos(angle)
    s = math.sin(angle)
    t = 1 - c
    x, y, z = axis
    
    return np.array([
        [t*x*x + c, t*x*y - s*z, t*x*z + s*y],
        [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
        [t*x*z - s*y, t*y*z + s*x, t*z*z + c]
    ])


def random_rotation_matrix():
    axis = [random.gauss(0, 1) for _ in range(3)]
    angle = random.uniform(0, 2 * math.pi)
    return rotation_matrix(axis, angle)


class DockingService:
    def __init__(self):
        self.scoring_weights = {
            'vdw': 0.4,
            'electrostatic': 0.3,
            'clash': 2.0,
            'hbond': -1.5,
            'hydrophobic': -0.5
        }
    
    def score_pose(self, ligand_atoms, receptor_atoms, binding_site):
        center = np.array(binding_site.center)
        radius = binding_site.radius
        
        score = 0.0
        clash_count = 0
        hbond_count = 0
        
        ligand_center = np.mean([[a['x'], a['y'], a['z']] for a in ligand_atoms], axis=0)
        dist_to_center = np.linalg.norm(ligand_center - center)
        
        if dist_to_center > radius:
            score += (dist_to_center - radius) * 5.0
        
        for latom in ligand_atoms:
            lpos = np.array([latom['x'], latom['y'], latom['z']])
            lcharge = ATOMIC_CHARGES.get(latom['element'], 0.0)
            lradius = VAN_DER_WAALS_RADII.get(latom['element'], 1.7)
            
            for ratom in receptor_atoms:
                rpos = np.array([ratom['x'], ratom['y'], ratom['z']])
                rcharge = ratom.get('charge', ATOMIC_CHARGES.get(ratom['element'], 0.0))
                rradius = VAN_DER_WAALS_RADII.get(ratom['element'], 1.7)
                
                dist = np.linalg.norm(lpos - rpos)
                
                if dist < 0.1:
                    dist = 0.1
                
                sigma = (lradius + rradius) / 2.0
                
                if dist < sigma * 0.8:
                    score += self.scoring_weights['clash'] * (sigma * 0.8 - dist)
                    clash_count += 1
                elif dist < sigma * 1.5:
                    ratio = (sigma / dist) ** 6
                    vdw_energy = 4.0 * ratio * (ratio - 1.0)
                    score += self.scoring_weights['vdw'] * vdw_energy
                
                if dist < 6.0 and lcharge != 0 and rcharge != 0:
                    electrostatic = lcharge * rcharge / (dist * 4.0)
                    score += self.scoring_weights['electrostatic'] * electrostatic
                
                hbond_donors = {'N', 'O'}
                if (latom['element'] in hbond_donors and ratom['element'] in hbond_donors 
                    and 2.5 < dist < 3.5):
                    angle_factor = 1.0
                    hbond_count += 1
                    score += self.scoring_weights['hbond'] * angle_factor
        
        return score, clash_count, hbond_count
    
    def apply_transformation(self, atoms, rot_matrix, translation):
        transformed = []
        for atom in atoms:
            pos = np.array([atom['x'], atom['y'], atom['z']])
            new_pos = rot_matrix @ pos + np.array(translation)
            new_atom = dict(atom)
            new_atom['x'] = float(new_pos[0])
            new_atom['y'] = float(new_pos[1])
            new_atom['z'] = float(new_pos[2])
            transformed.append(new_atom)
        return transformed
    
    def dock(self, ligand_atoms, receptor_atoms, binding_site,
             num_poses=10, num_rotations=60, num_translations=20,
             callback=None):
        center = np.array(binding_site.center)
        radius = binding_site.radius
        
        ligand_center = np.mean([[a['x'], a['y'], a['z']] for a in ligand_atoms], axis=0)
        
        poses = []
        total_iterations = num_rotations * num_translations
        iteration = 0
        
        for rot_idx in range(num_rotations):
            rot_mat = random_rotation_matrix()
            
            for trans_idx in range(num_translations):
                iteration += 1
                
                r = random.uniform(0, radius * 0.5)
                theta = random.uniform(0, 2 * math.pi)
                phi = random.uniform(0, math.pi)
                
                translation = center - ligand_center + np.array([
                    r * math.sin(phi) * math.cos(theta),
                    r * math.sin(phi) * math.sin(theta),
                    r * math.cos(phi)
                ])
                
                transformed = self.apply_transformation(ligand_atoms, rot_mat, translation)
                
                score, clash_count, hbond_count = self.score_pose(
                    transformed, receptor_atoms, binding_site
                )
                
                pose = DockingPose(
                    atoms=transformed,
                    score=score,
                    rotation=rot_mat.flatten().tolist(),
                    translation=translation.tolist(),
                    clash_count=clash_count,
                    hbond_count=hbond_count
                )
                poses.append(pose)
        
        poses.sort(key=lambda p: p.score)
        
        unique_poses = []
        for pose in poses:
            if len(unique_poses) >= num_poses:
                break
            is_duplicate = False
            for existing in unique_poses:
                pos_diff = np.linalg.norm(
                    np.array(pose.translation) - np.array(existing.translation)
                )
                if pos_diff < 1.0:
                    is_duplicate = True
                    break
            if not is_duplicate and pose.clash_count < len(ligand_atoms) * 0.5:
                unique_poses.append(pose)
        
        for i, pose in enumerate(unique_poses):
            pose.rank = i + 1
        
        return unique_poses
    
    def get_preset_receptor(self, preset_name):
        if preset_name not in PRESET_RECEPTORS:
            return None
        
        preset = PRESET_RECEPTORS[preset_name]
        if preset["generate"] and not preset["atoms"]:
            pocket_type = preset.get("pocket_type", "hydrophobic")
            center = tuple(preset["binding_site"]["center"])
            preset["atoms"] = generate_pocket_atoms(pocket_type, center)
        
        return preset


docking_service = DockingService()
