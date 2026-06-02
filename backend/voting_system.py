import json
from collections import defaultdict

class VotingSystem:
    def __init__(self):
        self.conformations = defaultdict(dict)
        self.votes = defaultdict(lambda: defaultdict(set))
        self.user_votes = defaultdict(dict)
    
    def add_conformation(self, session_id, conformation_id, conformation_data):
        self.conformations[session_id][conformation_id] = conformation_data
    
    def get_conformations(self, session_id):
        return list(self.conformations[session_id].values())
    
    def vote(self, session_id, conformation_id, user_id):
        if conformation_id not in self.conformations[session_id]:
            return False
        
        if user_id in self.user_votes[session_id]:
            old_conf_id = self.user_votes[session_id][user_id]
            if old_conf_id == conformation_id:
                self.votes[session_id][old_conf_id].discard(user_id)
                del self.user_votes[session_id][user_id]
                return True
            self.votes[session_id][old_conf_id].discard(user_id)
        
        self.user_votes[session_id][user_id] = conformation_id
        self.votes[session_id][conformation_id].add(user_id)
        return True
    
    def get_vote_count(self, session_id, conformation_id):
        return len(self.votes[session_id].get(conformation_id, set()))
    
    def get_all_votes(self, session_id):
        result = {}
        for conf_id, voters in self.votes[session_id].items():
            result[conf_id] = {
                'count': len(voters),
                'voters': list(voters)
            }
        return result
    
    def get_best_conformation(self, session_id):
        if not self.votes[session_id]:
            return None
        
        max_votes = -1
        best_conf_id = None
        
        for conf_id, voters in self.votes[session_id].items():
            if len(voters) > max_votes:
                max_votes = len(voters)
                best_conf_id = conf_id
        
        if best_conf_id:
            return {
                'conformation_id': best_conf_id,
                'vote_count': max_votes,
                'conformation': self.conformations[session_id].get(best_conf_id)
            }
        return None
    
    def get_user_vote(self, session_id, user_id):
        return self.user_votes[session_id].get(user_id)
    
    def clear_session(self, session_id):
        if session_id in self.conformations:
            del self.conformations[session_id]
        if session_id in self.votes:
            del self.votes[session_id]
        if session_id in self.user_votes:
            del self.user_votes[session_id]
