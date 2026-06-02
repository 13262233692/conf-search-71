from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_socketio import SocketIO, emit
from voting_system import VotingSystem
from docking_service import docking_service, BindingSite, generate_pocket_atoms, PRESET_RECEPTORS
import uuid
import json

try:
    from conformation_search import ConformationSearcher
    USE_REAL_RDKIT = True
except ImportError:
    from conformation_search_mock import MockConformationSearcher as ConformationSearcher
    USE_REAL_RDKIT = False
    print("Warning: RDKit not found, using mock conformation searcher")

app = Flask(__name__)
CORS(app)
socketio = SocketIO(app, cors_allowed_origins="*", async_mode='eventlet')

searcher = ConformationSearcher()
voting_system = VotingSystem()

current_sessions = {}

@app.route('/api/health', methods=['GET'])
def health_check():
    return jsonify({"status": "healthy", "service": "conformation-search"})

@app.route('/api/search', methods=['POST'])
def start_conformation_search():
    data = request.json
    smiles = data.get('smiles')
    method = data.get('method', 'etkdg')
    num_confs = data.get('num_confs', 10)
    session_id = data.get('session_id', str(uuid.uuid4()))
    
    if not smiles:
        return jsonify({"error": "SMILES is required"}), 400
    
    def search_callback(conformation, index, total):
        conf_id = f"{session_id}_{index}"
        conformation['id'] = conf_id
        conformation['session_id'] = session_id
        conformation['index'] = index
        
        socketio.emit('conformation_update', {
            'session_id': session_id,
            'conformation': conformation,
            'progress': (index + 1) / total * 100,
            'current': index + 1,
            'total': total
        })
        
        voting_system.add_conformation(session_id, conf_id, conformation)
    
    try:
        results = searcher.search_conformations(
            smiles=smiles,
            method=method,
            num_confs=num_confs,
            callback=search_callback
        )
        
        current_sessions[session_id] = {
            'smiles': smiles,
            'method': method,
            'results': results,
            'status': 'completed'
        }
        
        socketio.emit('search_complete', {
            'session_id': session_id,
            'total_conformations': len(results),
            'results': results
        })
        
        return jsonify({
            'session_id': session_id,
            'status': 'completed',
            'total_conformations': len(results)
        })
        
    except Exception as e:
        socketio.emit('search_error', {
            'session_id': session_id,
            'error': str(e)
        })
        return jsonify({"error": str(e)}), 500

@socketio.on('connect')
def handle_connect():
    print('Client connected')
    emit('connected', {'status': 'connected'})

@socketio.on('disconnect')
def handle_disconnect():
    print('Client disconnected')

@socketio.on('join_session')
def handle_join_session(data):
    session_id = data.get('session_id')
    if session_id:
        conformations = voting_system.get_conformations(session_id)
        votes = voting_system.get_all_votes(session_id)
        emit('session_data', {
            'session_id': session_id,
            'conformations': conformations,
            'votes': votes
        })

@socketio.on('vote')
def handle_vote(data):
    session_id = data.get('session_id')
    conformation_id = data.get('conformation_id')
    user_id = data.get('user_id', 'anonymous')
    
    success = voting_system.vote(session_id, conformation_id, user_id)
    
    if success:
        all_votes = voting_system.get_all_votes(session_id)
        best_conf = voting_system.get_best_conformation(session_id)
        
        socketio.emit('vote_update', {
            'session_id': session_id,
            'votes': all_votes,
            'best_conformation': best_conf
        })

@app.route('/api/votes/<session_id>', methods=['GET'])
def get_votes(session_id):
    votes = voting_system.get_all_votes(session_id)
    best = voting_system.get_best_conformation(session_id)
    return jsonify({
        'votes': votes,
        'best_conformation': best
    })

@app.route('/api/vote', methods=['POST'])
def post_vote():
    data = request.json
    session_id = data.get('session_id')
    conformation_id = data.get('conformation_id')
    user_id = data.get('user_id', 'anonymous')
    
    success = voting_system.vote(session_id, conformation_id, user_id)
    
    if success:
        socketio.emit('vote_update', {
            'session_id': session_id,
            'votes': voting_system.get_all_votes(session_id),
            'best_conformation': voting_system.get_best_conformation(session_id)
        })
        return jsonify({'success': True})
    return jsonify({'success': False}), 400

@app.route('/api/docking/presets', methods=['GET'])
def get_docking_presets():
    presets = []
    for key, preset in PRESET_RECEPTORS.items():
        presets.append({
            'id': key,
            'name': preset['name'],
            'description': preset['description'],
            'binding_site': preset['binding_site'],
            'pocket_type': preset.get('pocket_type', 'hydrophobic')
        })
    return jsonify(presets)

@app.route('/api/docking/preset/<preset_id>', methods=['GET'])
def get_docking_preset(preset_id):
    preset = docking_service.get_preset_receptor(preset_id)
    if preset is None:
        return jsonify({"error": "Preset not found"}), 404
    
    if not preset['atoms']:
        pocket_type = preset.get('pocket_type', 'hydrophobic')
        center = tuple(preset['binding_site']['center'])
        preset['atoms'] = generate_pocket_atoms(pocket_type, center)
    
    return jsonify({
        'id': preset_id,
        'name': preset['name'],
        'description': preset['description'],
        'binding_site': preset['binding_site'],
        'atoms': preset['atoms']
    })

@app.route('/api/docking/dock', methods=['POST'])
def run_docking():
    data = request.json
    
    ligand_atoms = data.get('ligand_atoms')
    receptor_preset = data.get('receptor_preset', 'hydrophobic_pocket')
    receptor_atoms_data = data.get('receptor_atoms')
    binding_site_data = data.get('binding_site', {'center': [0, 0, 0], 'radius': 8.0})
    num_poses = data.get('num_poses', 5)
    num_rotations = data.get('num_rotations', 60)
    num_translations = data.get('num_translations', 20)
    session_id = data.get('session_id', str(uuid.uuid4()))
    
    if not ligand_atoms:
        return jsonify({"error": "Ligand atoms are required"}), 400
    
    try:
        if receptor_atoms_data:
            receptor_atoms = receptor_atoms_data
        else:
            preset = docking_service.get_preset_receptor(receptor_preset)
            if preset is None:
                return jsonify({"error": f"Unknown preset: {receptor_preset}"}), 400
            if not preset['atoms']:
                pocket_type = preset.get('pocket_type', 'hydrophobic')
                center = tuple(preset['binding_site']['center'])
                preset['atoms'] = generate_pocket_atoms(pocket_type, center)
            receptor_atoms = preset['atoms']
        
        binding_site = BindingSite(
            center=binding_site_data.get('center', [0, 0, 0]),
            radius=binding_site_data.get('radius', 8.0),
            name=binding_site_data.get('name', 'default')
        )
        
        socketio.emit('docking_progress', {
            'session_id': session_id,
            'status': 'running',
            'message': 'Starting docking simulation...'
        })
        
        poses = docking_service.dock(
            ligand_atoms=ligand_atoms,
            receptor_atoms=receptor_atoms,
            binding_site=binding_site,
            num_poses=num_poses,
            num_rotations=num_rotations,
            num_translations=num_translations
        )
        
        results = []
        for pose in poses:
            results.append({
                'rank': pose.rank,
                'score': pose.score,
                'clash_count': pose.clash_count,
                'hbond_count': pose.hbond_count,
                'atoms': pose.atoms,
                'rotation': pose.rotation,
                'translation': pose.translation
            })
        
        socketio.emit('docking_complete', {
            'session_id': session_id,
            'poses': results,
            'receptor_atoms': receptor_atoms,
            'binding_site': {
                'center': binding_site.center,
                'radius': binding_site.radius
            }
        })
        
        return jsonify({
            'session_id': session_id,
            'poses': results,
            'total_poses': len(results)
        })
        
    except Exception as e:
        socketio.emit('docking_error', {
            'session_id': session_id,
            'error': str(e)
        })
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    socketio.run(app, host='0.0.0.0', port=5000, debug=True)
