const { createApp } = Vue;

createApp({
    data() {
        return {
            connected: false,
            socket: null,
            currentSessionId: null,
            joinSessionId: '',
            userId: 'user_' + Math.random().toString(36).substr(2, 9),
            
            smilesInput: '',
            searchMethod: 'etkdg',
            numConformations: 10,
            
            isSearching: false,
            searchProgress: 0,
            currentConformation: 0,
            totalConformations: 0,
            
            conformations: [],
            selectedConformation: null,
            votes: {},
            userVote: null,
            bestConformation: null,
            detectedRingInfo: null,
            
            isDocking: false,
            dockingReceptor: 'hydrophobic_pocket',
            dockingRadius: 8.0,
            dockingNumPoses: 5,
            dockingRotations: 60,
            dockingResult: null,
            selectedDockingPose: null,
            showDockingView: false,
            
            viewer: null,
            spinning: false,
            currentStyle: 'stick'
        };
    },
    
    computed: {
        sortedConformations() {
            return [...this.conformations].sort((a, b) => {
                const votesA = this.getVoteCount(a.id);
                const votesB = this.getVoteCount(b.id);
                if (votesB !== votesA) return votesB - votesA;
                return (a.rank || 0) - (b.rank || 0);
            });
        }
    },
    
    mounted() {
        this.initViewer();
        this.initSocket();
    },
    
    methods: {
        initViewer() {
            const config = { backgroundColor: '#000000' };
            this.viewer = $3Dmol.createViewer('molviewer', config);
            this.viewer.zoomTo();
            this.viewer.render();
        },
        
        initSocket() {
            this.socket = io('http://localhost:5000', {
                transports: ['websocket', 'polling']
            });
            
            this.socket.on('connect', () => {
                this.connected = true;
            });
            
            this.socket.on('disconnect', () => {
                this.connected = false;
            });
            
            this.socket.on('conformation_update', (data) => {
                this.searchProgress = data.progress;
                this.currentConformation = data.current;
                this.totalConformations = data.total;
                this.conformations.push(data.conformation);
                
                if (this.conformations.length === 1) {
                    this.selectConformation(data.conformation);
                }
            });
            
            this.socket.on('search_complete', (data) => {
                this.isSearching = false;
                this.searchProgress = 100;
            });
            
            this.socket.on('search_error', (data) => {
                this.isSearching = false;
                alert('搜索错误: ' + data.error);
            });
            
            this.socket.on('session_data', (data) => {
                this.conformations = data.conformations || [];
                this.votes = data.votes || {};
                if (this.conformations.length > 0 && !this.selectedConformation) {
                    this.selectConformation(this.conformations[0]);
                }
            });
            
            this.socket.on('vote_update', (data) => {
                this.votes = data.votes || {};
                this.bestConformation = data.best_conformation;
                this.updateUserVote();
            });
            
            this.socket.on('docking_progress', (data) => {
                this.isDocking = true;
            });
            
            this.socket.on('docking_complete', (data) => {
                this.isDocking = false;
                this.dockingResult = data;
                if (data.poses && data.poses.length > 0) {
                    this.selectedDockingPose = data.poses[0];
                    this.showDockingView = true;
                    this.renderDockingView(data.poses[0], data.receptor_atoms, data.binding_site);
                }
            });
            
            this.socket.on('docking_error', (data) => {
                this.isDocking = false;
                alert('对接错误: ' + data.error);
            });
        },
        
        setSmiles(smiles) {
            this.smilesInput = smiles;
            this.checkForLargeRing();
        },
        
        checkForLargeRing() {
            const smiles = this.smilesInput;
            if (!smiles) {
                this.detectedRingInfo = null;
                return;
            }
            
            const ringSizes = this.detectRingSizes(smiles);
            const largeRings = ringSizes.filter(s => s >= 10);
            
            if (largeRings.length > 0) {
                this.detectedRingInfo = {
                    num_large_rings: largeRings.length,
                    max_ring_size: Math.max(...ringSizes)
                };
            } else {
                this.detectedRingInfo = null;
            }
        },
        
        detectRingSizes(smiles) {
            const ringSizes = [];
            const ringMap = {};
            
            for (let i = 0; i < smiles.length; i++) {
                const char = smiles[i];
                
                if (char === '%' && i + 2 < smiles.length) {
                    const ringNum = parseInt(smiles.substring(i + 1, i + 3));
                    i += 2;
                    if (ringMap[ringNum] !== undefined) {
                        ringSizes.push(i - ringMap[ringNum] + 1);
                        delete ringMap[ringNum];
                    } else {
                        ringMap[ringNum] = i;
                    }
                } else if (char >= '0' && char <= '9') {
                    const ringNum = parseInt(char);
                    if (ringMap[ringNum] !== undefined) {
                        ringSizes.push(i - ringMap[ringNum] + 1);
                        delete ringMap[ringNum];
                    } else {
                        ringMap[ringNum] = i;
                    }
                }
            }
            
            return ringSizes;
        },
        
        async startSearch() {
            if (!this.smilesInput) return;
            
            this.isSearching = true;
            this.searchProgress = 0;
            this.conformations = [];
            this.votes = {};
            this.userVote = null;
            this.bestConformation = null;
            this.dockingResult = null;
            this.selectedDockingPose = null;
            this.showDockingView = false;
            this.currentSessionId = 'session_' + Date.now();
            
            try {
                const response = await fetch('http://localhost:5000/api/search', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        smiles: this.smilesInput,
                        method: this.searchMethod,
                        num_confs: parseInt(this.numConformations),
                        session_id: this.currentSessionId
                    })
                });
                
                const result = await response.json();
                if (!response.ok) {
                    throw new Error(result.error || 'Search failed');
                }
            } catch (error) {
                this.isSearching = false;
                alert('启动搜索失败: ' + error.message);
            }
        },
        
        selectConformation(conf) {
            this.selectedConformation = conf;
            if (!this.showDockingView) {
                this.renderMolecule(conf);
            }
            this.updateUserVote();
        },
        
        renderMolecule(conf) {
            if (!this.viewer || !conf) return;
            
            this.viewer.clear();
            
            let xyzData = '';
            if (conf.atoms) {
                xyzData = conf.atoms.length + '\n\n';
                conf.atoms.forEach(atom => {
                    xyzData += `${atom.element} ${atom.x} ${atom.y} ${atom.z}\n`;
                });
            }
            
            this.viewer.addModel(xyzData, 'xyz');
            this.applyStyle(this.currentStyle);
            this.viewer.zoomTo();
            this.viewer.render();
        },
        
        renderDockingView(pose, receptorAtoms, bindingSite) {
            if (!this.viewer || !pose) return;
            
            this.viewer.clear();
            
            let receptorXyz = receptorAtoms.length + '\n\n';
            receptorAtoms.forEach(atom => {
                receptorXyz += `${atom.element} ${atom.x} ${atom.y} ${atom.z}\n`;
            });
            
            const receptorModel = this.viewer.addModel(receptorXyz, 'xyz');
            this.viewer.setStyle({ model: receptorModel }, { 
                stick: { radius: 0.15, colorscheme: 'default' },
                sphere: { scale: 0.2, colorscheme: 'default', opacity: 0.6 }
            });
            
            if (bindingSite) {
                this.viewer.addSphere({
                    center: { x: bindingSite.center[0], y: bindingSite.center[1], z: bindingSite.center[2] },
                    radius: bindingSite.radius,
                    color: '#4488ff',
                    opacity: 0.08,
                    wireframe: true
                });
            }
            
            let ligandXyz = pose.atoms.length + '\n\n';
            pose.atoms.forEach(atom => {
                ligandXyz += `${atom.element} ${atom.x} ${atom.y} ${atom.z}\n`;
            });
            
            const ligandModel = this.viewer.addModel(ligandXyz, 'xyz');
            this.viewer.setStyle({ model: ligandModel }, { 
                stick: { radius: 0.25, color: '#00ff88' },
                sphere: { scale: 0.35, color: '#00ff88' }
            });
            
            this.viewer.zoomTo({ model: ligandModel });
            this.viewer.render();
        },
        
        toggleDockingView() {
            if (this.showDockingView) {
                this.showDockingView = false;
                if (this.selectedConformation) {
                    this.renderMolecule(this.selectedConformation);
                }
            } else {
                this.showDockingView = true;
                if (this.dockingResult && this.selectedDockingPose) {
                    this.renderDockingView(
                        this.selectedDockingPose, 
                        this.dockingResult.receptor_atoms,
                        this.dockingResult.binding_site
                    );
                }
            }
        },
        
        selectDockingPose(pose) {
            this.selectedDockingPose = pose;
            if (this.showDockingView && this.dockingResult) {
                this.renderDockingView(pose, this.dockingResult.receptor_atoms, this.dockingResult.binding_site);
            }
        },
        
        async startDocking() {
            if (!this.selectedConformation || !this.selectedConformation.atoms) return;
            
            this.isDocking = true;
            this.dockingResult = null;
            this.selectedDockingPose = null;
            
            try {
                const response = await fetch('http://localhost:5000/api/docking/dock', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        ligand_atoms: this.selectedConformation.atoms,
                        receptor_preset: this.dockingReceptor,
                        binding_site: {
                            center: [0, 0, 0],
                            radius: parseFloat(this.dockingRadius)
                        },
                        num_poses: parseInt(this.dockingNumPoses),
                        num_rotations: parseInt(this.dockingRotations),
                        num_translations: 20,
                        session_id: this.currentSessionId || 'session_' + Date.now()
                    })
                });
                
                const result = await response.json();
                if (!response.ok) {
                    throw new Error(result.error || 'Docking failed');
                }
            } catch (error) {
                this.isDocking = false;
                alert('对接失败: ' + error.message);
            }
        },
        
        clearDocking() {
            this.dockingResult = null;
            this.selectedDockingPose = null;
            this.showDockingView = false;
            if (this.selectedConformation) {
                this.renderMolecule(this.selectedConformation);
            }
        },
        
        setStyle(style) {
            this.currentStyle = style;
            if (this.showDockingView && this.dockingResult && this.selectedDockingPose) {
                this.renderDockingView(
                    this.selectedDockingPose, 
                    this.dockingResult.receptor_atoms,
                    this.dockingResult.binding_site
                );
                return;
            }
            if (!this.viewer) return;
            this.applyStyle(style);
            this.viewer.render();
        },
        
        applyStyle(style) {
            if (!this.viewer) return;
            this.viewer.setStyle({}, {});
            
            if (style === 'stick') {
                this.viewer.setStyle({}, { stick: { radius: 0.2, colorscheme: 'default' } });
            } else if (style === 'sphere') {
                this.viewer.setStyle({}, { sphere: { scale: 0.3, colorscheme: 'default' } });
            } else if (style === 'cartoon') {
                this.viewer.setStyle({}, { stick: { radius: 0.2 }, sphere: { scale: 0.2 } });
            }
        },
        
        spinToggle() {
            this.spinning = !this.spinning;
            if (this.spinning) {
                this.viewer.spin('y', 1);
            } else {
                this.viewer.spin(false);
            }
        },
        
        vote(conformationId) {
            this.socket.emit('vote', {
                session_id: this.currentSessionId,
                conformation_id: conformationId,
                user_id: this.userId
            });
        },
        
        getVoteCount(conformationId) {
            return this.votes[conformationId]?.count || 0;
        },
        
        updateUserVote() {
            for (const [confId, voteData] of Object.entries(this.votes)) {
                if (voteData.voters?.includes(this.userId)) {
                    this.userVote = confId;
                    return;
                }
            }
            this.userVote = null;
        },
        
        joinSession() {
            if (!this.joinSessionId) return;
            
            this.currentSessionId = this.joinSessionId;
            this.conformations = [];
            this.votes = {};
            
            this.socket.emit('join_session', {
                session_id: this.joinSessionId
            });
            
            this.joinSessionId = '';
        }
    }
}).mount('#app');
