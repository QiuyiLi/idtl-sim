from shared_utility import *
from species_tree import *

import os


class GeneTree(GenericTree):
    # static properties
    dup_recombination = 0
    hemiplasy = 0
    lambda_dup = None
    lambda_loss = None
    lambda_trans = None
    full_events = []

    def __init__(self,
                 time_sequences,
                 species_tree,
                 coalescent_process):
        GenericTree.__init__(self)
        self.time_sequences = time_sequences
        self.coalescent_process = coalescent_process
        self.construct_gene_nodes()

        max_node_id = -1
        for node in self.nodes:
            if (node.node_id > max_node_id):
                max_node_id = node.node_id
                self.root = node
            self.nodes_id_dict[node.node_id] = node
            self.nodes_name_dict[node.name] = node

        for node in self.nodes:
            clade = node.name.split('*')[:-1]
            clade = [int(j) for j in clade]
            node.clade = clade
            if (node.children and not node.clade_split):
                for i in range(len(node.children)):
                    node_name = self.node_by_id(node.children[i]).name 
                    clade = node_name.split('*')[:-1]
                    clade = [int(j) for j in clade]
                    node.clade_split.append(clade)
    
        self.species_tree = species_tree
        self.leaves = [node.node_id for node in self.nodes if not node.children]
        self.total_distance = self.distance_to_root_recurse(node_id=self.leaves[0])
        return

    # string replacement modified for the "*" representation
    def star_replace(self, string, substring):
        a = string.split('*')[:-1]
        b = substring.split('*')[:-1]
        diff = set(a).difference(set(b))
        return ''.join([e + '*' for e in sorted(list(diff))])

    def distance_from_to(self, node_name, parent_name):
        for leaf, sequence in self.time_sequences.items():
            if (node_name.count('*') == 1 and node_name[0] == leaf):
                for pair in sequence:
                    if (pair[0] == parent_name):
                        return pair[1]
            else:
                prev_pair = None
                for pair in sequence:
                    if (prev_pair != None and prev_pair[0] == node_name and pair[0] == parent_name):
                        return pair[1] - prev_pair[1]
                    prev_pair = pair
        return None

    def create_skbio_tree_recurse(self, skbio_tree_node):
        # one node (leaf)
        if (skbio_tree_node.name.count('*') == 1):
            skbio_tree_node.length = self.distance_from_to(skbio_tree_node.name, skbio_tree_node.parent.name)
            return
        # two nodes
        elif (len(skbio_tree_node.name) == 4):
            child_one_name = skbio_tree_node.name[:2]
            child_two_name = skbio_tree_node.name[2:]
            child_one = skbio.tree.TreeNode(name=child_one_name, 
                                            length=self.distance_from_to(child_one_name, skbio_tree_node.name), 
                                            parent=skbio_tree_node)
            child_two = skbio.tree.TreeNode(name=child_two_name, 
                                            length=self.distance_from_to(child_two_name, skbio_tree_node.name),
                                            parent=skbio_tree_node)
            skbio_tree_node.children = [child_one, child_two]
            return
        is_found = False
        for leaf, sequence in self.time_sequences.items():
            prev_pair = None
            for pair in sequence:
                if (prev_pair != None and skbio_tree_node.name == pair[0]):
                    child_one_name = prev_pair[0]
                    child_two_name = self.star_replace(skbio_tree_node.name, prev_pair[0])
                    child_one = skbio.tree.TreeNode(name=child_one_name, 
                                                    length=self.distance_from_to(child_one_name, skbio_tree_node.name), 
                                                    parent=skbio_tree_node)
                    child_two = skbio.tree.TreeNode(name=child_two_name, 
                                                    length=self.distance_from_to(child_two_name, skbio_tree_node.name),
                                                    parent=skbio_tree_node)
                    self.create_skbio_tree_recurse(child_one)
                    self.create_skbio_tree_recurse(child_two)
                    skbio_tree_node.children = [child_one, child_two]
                    is_found = True
                    break
                prev_pair = pair
            if (is_found):
                break
        return

    # construct ghe gene tree in newick format from time sequence
    def construct_gene_nodes(self):
        tree = skbio.tree.TreeNode() # empty tree, initialization
        time_seq = self.time_sequences.copy()
        for k, v in self.time_sequences.items():
            if (not v): del time_seq[k]
        if (len(time_seq) > 0 and next(iter(time_seq.values()))):
            tree.name = next(iter(time_seq.values()))[-1][0]
            self.create_skbio_tree_recurse(tree)
            tree.length = None
            self.skbio_tree = tree
            output_path = 'output/temp_gene_nodes_table.txt'
            super().newick_to_table(skbio_tree=tree, output_path=output_path)
            super().construct_nodes(path=output_path, process_tree=False)   
        else: 
            tree.name = next(iter(self.time_sequences)) + '*' # empty
            tree.length = None
            self.skbio_tree = tree

            self.nodes.append(TreeNode(node_id=0,
                                       name=tree.name,
                                       parent=-1,
                                       distance_to_parent=-1.0))
        return

    def construct_final_gene_nodes(self, tree):
        self.nodes = []
        self.nodes_id_dict = {}
        self.nodes_name_dict = {}
        output_path = 'output/temp_final_gene_nodes_table.txt'
        self.skbio_tree = super().newick_to_table(skbio_tree=tree, output_path=output_path)
        super().construct_nodes(path=output_path, process_tree=True)
        max_node_id = -1
        for node in self.nodes:
            if (node.node_id > max_node_id):
                max_node_id = node.node_id
                self.root = node
            self.nodes_id_dict[node.node_id] = node
            self.nodes_name_dict[node.name] = node
        return

    def get_lambda_dup(self, clade):
            indices = []
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
            return mean(GeneTree.lambda_dup[indices])

    def get_lambda_loss(self, clade):
            indices = []
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
            return mean(GeneTree.lambda_loss[indices])

    def get_lambda_trans(self, clade):
            indices = []
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
            return mean(GeneTree.lambda_trans[indices])

    def map_gene_id_to_species_id(self, gene_id):
        coalescent_process = self.coalescent_process
        species_id = None
        gene_name = self.node_by_id(gene_id).name
        if (coalescent_process):        # non-trivial case
            for k, v in coalescent_process.items():
                for elem in v:
                    if (gene_name in elem['to_set'] and gene_name not in elem['from_set']):
                        species_id = int(k)
            if (species_id == None):
                species_id = int(gene_name[:-1])
        else:                           # trivial case
            species_id = int(gene_name[:-1])     
        return species_id

    def map_event_to_species_tree(self, gene_id, event_height, species_id=None):
        if (species_id == None):
            species_id = self.map_gene_id_to_species_id(gene_id)
        distance_above_species_node = event_height - SpeciesTree.global_species_tree.distance_to_leaf(species_id, 0)
        if (species_id == SpeciesTree.global_species_tree.root.node_id):
            return species_id, distance_above_species_node
        else:
            species_id_parent = SpeciesTree.global_species_tree.node_by_id(species_id).parent_id
            species_distance_parent = SpeciesTree.global_species_tree.distance_to_leaf(species_id_parent, 0)
            if (species_distance_parent > event_height):
                return species_id, distance_above_species_node
            else:
                return self.map_event_to_species_tree(gene_id, event_height, species_id=species_id_parent)

    def find_trans_target(self, event_height, node_id):
        tree = SpeciesTree.global_species_tree
        species_nodes = tree.nodes
        origin_species_id = self.map_gene_id_to_species_id(node_id)
        nodes_list = []
        for node in species_nodes:
            if (node.node_id == origin_species_id):
                continue
            if (node.node_id == tree.root.node_id):
                continue
            parent_height = tree.distance_to_leaf(tree.node_by_id(node.parent_id).node_id, 0)
            if (parent_height > event_height):
                node_height = tree.distance_to_leaf(node.node_id, 0)
                if (node_height <= event_height):
                    nodes_list.append(node.node_id)
        return np.random.choice(nodes_list), origin_species_id

    # # find the points of duplicatons and losses recursively
    # def dlt_process_recurse(self, tree, distance, events):
    #     node = self.nodes_name_dict[tree.name]
    #     distance_dup = np.random.exponential(scale=1.0/self.get_lambda_dup(node.name))
    #     distance_loss = np.random.exponential(scale=1.0/self.get_lambda_loss(node.name))
    #     distance_trans = np.random.exponential(scale=1.0/self.get_lambda_trans(node.name))
    #     if (distance_dup < min(distance_loss, distance_trans) and distance_dup < distance):      # duplication happens first
    #         Debug.log(header='duplication at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                   ' with distance ' + str(distance - distance_dup) + '\n')
    #         event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_dup
    #         species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                event_height=event_height, 
    #                                                                                species_id=None)
    #         events.append({
    #             'type': 'duplication',
    #             'gene_node_id': node.node_id, 
    #             'gene_node_name': node.name, 
    #             'distance_to_gene_node': distance - distance_dup,
    #             'event_height': event_height,
    #             'species_node_id': species_id,
    #             'distance_to_species_node': distance_above_species_node,
    #             'index': -1
    #         })
    #         # looking for more events on the same branch
    #         self.dt_process_recurse(tree, distance - distance_dup, events) 
    #     elif (distance_trans <= min(distance_dup, distance_loss) and distance_trans < distance):
    #         event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_trans
    #         species_tree_height = SpeciesTree.global_species_tree.total_distance
    #         if (event_height < species_tree_height):
    #             Debug.log(header='transfer at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                       ' with distance ' + str(distance - distance_trans) + '\n')
    #             target, origin_species_id = self.find_trans_target(event_height, node.node_id)
    #             species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                    event_height=event_height, 
    #                                                                                    species_id=None)
    #             if(target):
    #                 events.append({
    #                     'type': 'transfer',
    #                     'gene_node_id': node.node_id, 
    #                     'gene_node_name': node.name, 
    #                     'distance_to_gene_node': distance - distance_trans,
    #                     'target_species_id': target,
    #                     'event_height': event_height,
    #                     'species_node_id': species_id,
    #                     'distance_to_species_node': distance_above_species_node,
    #                     'index': -1
    #                 })
    #         self.dlt_process_recurse(tree, distance - distance_trans, events)
    #     elif (distance_loss <= min(distance_dup, distance_trans) and distance_loss < distance):      
    #         # loss happens first, the seaching process stops at the loss point
    #         event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_loss
    #         Debug.log(header='loss at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                   ' with distance ' + str(distance - distance_loss) + '\n')
    #         species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                event_height=event_height, 
    #                                                                                species_id=None)
    #         # species_id = self.node_by_id(species_id).parent_id
    #         events.append({
    #             'type': 'loss',
    #             'gene_node_id': node.node_id, 
    #             'gene_node_name': node.name, 
    #             'distance_to_gene_node': distance - distance_loss,
    #             'event_height': event_height,
    #             'species_node_id': species_id,
    #             'distance_to_species_node': distance_above_species_node,
    #             'index': -1
    #         })
    #     else:   # reach the end the current branch, looking for events in the 2 children branches
    #         Debug.log(header='nothing happened at node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
    #         if (node.children):     # if children branches exist
    #             child_one = tree.children[0]
    #             child_two = tree.children[1]
    #             distance_to_child_one = node.distance_to_children[0]
    #             distance_to_child_two = node.distance_to_children[1]
    #             self.dlt_process_recurse(child_one, distance_to_child_one, events)
    #             self.dlt_process_recurse(child_two, distance_to_child_two, events)
    #         else:       # if not exist, reach the leaves of the tree, searching process stops
    #             Debug.log(header='reach the end of node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
    #     return

    # # find the points of duplicatons and losses recursively
    # def dt_process_recurse(self, tree, distance, events):
    #     node = self.nodes_name_dict[tree.name]
    #     distance_dup = np.random.exponential(scale=1.0/self.get_lambda_dup(node.name))
    #     distance_loss = 10000
    #     distance_trans = np.random.exponential(scale=1.0/self.get_lambda_trans(node.name))
    #     if (distance_dup < min(distance_loss, distance_trans) and distance_dup < distance):      # duplication happens first
    #         Debug.log(header='duplication at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                   ' with distance ' + str(distance - distance_dup) + '\n')
    #         event_height = event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_dup
    #         species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                event_height=event_height, 
    #                                                                                species_id=None)
    #         events.append({
    #             'type': 'duplication',
    #             'gene_node_id': node.node_id, 
    #             'gene_node_name': node.name, 
    #             'distance_to_gene_node': distance - distance_dup,
    #             'event_height': event_height,
    #             'species_node_id': species_id,
    #             'distance_to_species_node': distance_above_species_node,
    #             'index': -1
    #         })
    #         self.dt_process_recurse(tree, distance - distance_dup, events) # looking for more events on the same branch
    #     elif (distance_trans <= min(distance_dup, distance_loss) and distance_trans < distance):
    #         event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_trans
    #         species_tree_height = SpeciesTree.global_species_tree.total_distance
    #         if (event_height < species_tree_height):
    #             Debug.log(header='transfer at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                       ' with distance ' + str(distance - distance_trans) + '\n')
    #             target, origin_species_id = self.find_trans_target(event_height, node.node_id)
    #             species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                    event_height=event_height, 
    #                                                                                    species_id=None)
    #             if (target):
    #                 events.append({
    #                     'type': 'transfer',
    #                     'gene_node_id': node.node_id, 
    #                     'gene_node_name': node.name, 
    #                     'distance_to_gene_node': distance - distance_trans,
    #                     'target_species_id': target,
    #                     'event_height': event_height,
    #                     'species_node_id': species_id,
    #                     'distance_to_species_node': distance_above_species_node,
    #                     'index': -1
    #                 })
    #         self.dlt_process_recurse(tree, distance - distance_trans, events)
    #     elif (distance_loss <= min(distance_dup, distance_trans) and distance_loss < distance): 
    #         # loss happens first, the seaching process stops at the loss point
    #         Debug.log(header='loss at node ' + str(node.node_id) + ' (' + node.name + ')' + 
    #                   ' with distance ' + str(distance - distance_loss) + '\n')
    #         species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
    #                                                                                event_height=event_height, 
    #                                                                                species_id=None)
    #         # species_id = self.node_by_id(species_id).parent_id
    #         events.append({
    #             'type': 'loss',
    #             'gene_node_id': node.node_id, 
    #             'gene_node_name': node.name, 
    #             'distance_to_gene_node': distance - distance_loss,
    #             'event_height': event_height,
    #             'species_node_id': species_id, 
    #             'distance_to_species_node': distance_above_species_node,
    #             'index': -1
    #         })
    #     else:   # reach the end the current branch, looking for events in the 2 children branches
    #         Debug.log(header='nothing happened at node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
    #         if (node.children):     # if children branches exist
    #             child_one = tree.children[0]
    #             child_two = tree.children[1]
    #             distance_to_child_one = node.distance_to_children[0]
    #             distance_to_child_two = node.distance_to_children[1]
    #             self.dlt_process_recurse(child_one, distance_to_child_one, events)
    #             self.dlt_process_recurse(child_two, distance_to_child_two, events)
    #         else:       # if not exist, reach the leaves of the tree, searching process stops
    #             Debug.log(header='reach the end of node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
    #     return
    
    # # store the duplication events
    # def dlt_process(self, distance, event=None):
    #     events = []

    #     if (len(self.nodes) == 1):
    #         distance = event['distance_to_gene_node']
    #     self.dt_process_recurse(self.skbio_tree, 
    #                                 distance=distance, 
    #                                 events=events)  
    #     return events


       # find the points of duplicatons and losses recursively
    def dlt_process_recurse(self, tree, distance, events):
        node = self.nodes_name_dict[tree.name]
        distance_dup = np.random.exponential(scale=1.0/self.get_lambda_dup(node.name))
        distance_loss = np.random.exponential(scale=1.0/self.get_lambda_loss(node.name))
        distance_trans = np.random.exponential(scale=1.0/self.get_lambda_trans(node.name))
        if (distance_dup < min(distance_loss, distance_trans) and distance_dup < distance):      # duplication happens first
            Debug.log(header='duplication at node ' + str(node.node_id) + ' (' + node.name + ')' + 
                      ' with distance ' + str(distance - distance_dup) + '\n')
            event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_dup
            species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
                                                                                   event_height=event_height, 
                                                                                   species_id=None)
            events.append({
                'type': 'duplication',
                'gene_node_id': node.node_id, 
                'gene_node_name': node.name, 
                'distance_to_gene_node': distance - distance_dup,
                'event_height': event_height,
                'species_node_id': species_id,
                'distance_to_species_node': distance_above_species_node,
                'index': -1
            })
            # looking for more events on the same branch
            self.dlt_process_recurse(tree, distance - distance_dup, events) 
        elif (distance_trans <= min(distance_dup, distance_loss) and distance_trans < distance):
            event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_trans
            species_tree_height = SpeciesTree.global_species_tree.total_distance
            if (event_height < species_tree_height):
                Debug.log(header='transfer at node ' + str(node.node_id) + ' (' + node.name + ')' + 
                          ' with distance ' + str(distance - distance_trans) + '\n')
                target, origin_species_id = self.find_trans_target(event_height, node.node_id)
                species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
                                                                                       event_height=event_height, 
                                                                                       species_id=None)
                if(target):
                    events.append({
                        'type': 'transfer',
                        'gene_node_id': node.node_id, 
                        'gene_node_name': node.name, 
                        'distance_to_gene_node': distance - distance_trans,
                        'target_species_id': target,
                        'event_height': event_height,
                        'species_node_id': species_id,
                        'distance_to_species_node': distance_above_species_node,
                        'index': -1
                    })
            self.dlt_process_recurse(tree, distance - distance_trans, events)
        elif (distance_loss <= min(distance_dup, distance_trans) and distance_loss < distance):      
            # loss happens first, the seaching process stops at the loss point
            event_height = super().distance_to_leaf(node.node_id, 0) + distance - distance_loss
            Debug.log(header='loss at node ' + str(node.node_id) + ' (' + node.name + ')' + 
                      ' with distance ' + str(distance - distance_loss) + '\n')
            species_id, distance_above_species_node = self.map_event_to_species_tree(gene_id=node.node_id, 
                                                                                   event_height=event_height, 
                                                                                   species_id=None)
            # species_id = self.node_by_id(species_id).parent_id
            events.append({
                'type': 'loss',
                'gene_node_id': node.node_id, 
                'gene_node_name': node.name, 
                'distance_to_gene_node': distance - distance_loss,
                'event_height': event_height,
                'species_node_id': species_id,
                'distance_to_species_node': distance_above_species_node,
                'index': -1
            })
        else:   # reach the end the current branch, looking for events in the 2 children branches
            Debug.log(header='nothing happened at node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
            if (node.children):     # if children branches exist
                child_one = tree.children[0]
                child_two = tree.children[1]
                distance_to_child_one = node.distance_to_children[0]
                distance_to_child_two = node.distance_to_children[1]
                self.dlt_process_recurse(child_one, distance_to_child_one, events)
                self.dlt_process_recurse(child_two, distance_to_child_two, events)
            else:       # if not exist, reach the leaves of the tree, searching process stops
                Debug.log(header='reach the end of node ' + str(node.node_id) + ' (' + node.name + ')' + '\n')
        return
    
    # store the duplication events
    def dlt_process(self, distance, event=None):
        events = []

        if (len(self.nodes) == 1):
            distance = event['distance_to_gene_node']
        self.dlt_process_recurse(self.skbio_tree, 
                                    distance=distance, 
                                    events=events)  
        return events

    def find_ils(self, path):
        for i in range(len(self.nodes)):
            j = len(self.nodes)-1-i
            gene_node = self.node_by_id(j)
            gene_clade = gene_node.clade
            gene_splits = gene_node.clade_split
            find_species_node = False
            # find lowest possibile species node, then trace back
            for node in self.species_tree.nodes: 
                if (set(node.clade).issuperset(set(gene_clade))):
                    find_species_node = True
                    species_node = node
                    break
            if not find_species_node:
                continue
            species_clade = species_node.clade
            species_splits = species_node.clade_split
            find_ils = False
            if (gene_splits):
                gene_split_0 = set(gene_splits[0])
                gene_split_1 = set(gene_splits[1])
                for species_split in species_splits:
                    if (set(species_split).intersection(gene_split_0) and set(species_split).intersection(gene_split_1)):
                        find_ils = True
                        break
            if (find_ils):
                species_id = self.map_gene_id_to_species_id(j)
                index = Utility.increment()
                self.full_events.append({
                    'type': 'ils',
                    'gene_node_id': j, 
                    'gene_node_name': node.name, 
                    'species_node_id': species_id,
                    'index': index
                })
                # Debug.event_count['I'] += 1
                file_name = 'ils_' + str(index)
                f = open(os.path.join(path, file_name), 'w')
                f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
                f.close()
                # print('find ils at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
                #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))
            else: 
                if (gene_splits):
                    species_id = self.map_gene_id_to_species_id(j)
                    index = Utility.increment()
                    self.full_events.append({
                        'type': 's',
                        'gene_node_id': j, 
                        'gene_node_name': node.name, 
                        'species_node_id': species_id,
                        'index': index
                    })
                    file_name = 's_' + str(index)
                    f = open(os.path.join(path, file_name), 'w')
                    f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
                    f.close()
                    # print('find speciation at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
                    #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))

    # find the duplication subtree and do subtree coalescence
    def dt_subtree_recurse(self, event, node_id, coal_distance, path):

        if (event['type'] == 'transfer'): # node_id = target_id
            Debug.log(header='\n\n\n' + '='*80 + '\nCurrent event:' + '\n',
                      bodies=[event], pformat=True)
            tree = SpeciesTree.global_species_tree
            species_skbio_tree = tree.skbio_tree
            name = tree.nodes_id_dict[node_id].name

            subtree = species_skbio_tree.find(name).deepcopy()
            subtree_names = [node.name for node in subtree.traverse()]
            subtree_nodes = [node for node in tree.nodes if node.name in subtree_names]
            species_subtree = SpeciesTree(nodes=subtree_nodes)
            species_subtree.skbio_tree = subtree
            Debug.log(header='\nspecies_subtree_nodes:\n', bodies=species_subtree.nodes)

            distance_above_root = coal_distance
            Debug.log(header='\nspecies_subtree_coal:\n')
            if (GeneTree.hemiplasy == 1):
                species_subtree_coal_process, chosen_gene_name = species_subtree.incomplete_coalescent(distance_above_root)
            elif (GeneTree.hemiplasy == 0):
                species_subtree_coal_process = species_subtree.bounded_coalescent(distance_above_root)

            Debug.log(header='\nspecies_subtree_coal_process:\n', 
                      bodies=[species_subtree_coal_process], pformat=True)

            species_subtree_time_seq = species_subtree.time_sequences(coalescent_process=species_subtree_coal_process)
            Debug.log(header='\nspecies_subtree_time_seq:' + '\n',
                      bodies=[species_subtree_time_seq], pformat=True)

            # debug
            # Debug.save_tree_nodes(nodes=species_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'trans', node_id, distance_above_root), 
            #                       distance=distance_above_root)
            
            gene_subtree = GeneTree(time_sequences=species_subtree_time_seq, 
                                    species_tree=species_subtree, 
                                    coalescent_process=species_subtree_coal_process)
            gene_subtree.skbio_tree.length = event['event_height'] - gene_subtree.total_distance
            Debug.log(header='\ngene_subtree nodes:\n', bodies=gene_subtree.nodes)
            # debug
            # Debug.save_tree_nodes(nodes=gene_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'trans', node_id, distance_above_root), 
            #                       mode='a')

            # Debug.save_output(contents=[gene_subtree.skbio_tree],
            #                   path=Debug.subtree_file_name('output/newick_gene_subtrees', 'trans', node_id, distance_above_root))

            Debug.log(header='\ngene_subtree dlt_process:\n')
            gene_subtree_height = gene_subtree.total_distance
            gene_subtree_events = gene_subtree.dlt_process(event=event, distance=event['event_height'] - gene_subtree_height)
            Debug.log(header='\ngene_subtree events:\n', bodies=[gene_subtree_events], pformat=True)
            
            index = Utility.increment()
            event['index'] = index
            _id = 'trans_subtree_' + str(index)
            next_dir = os.path.join(path, _id)
            os.mkdir(next_dir)
            file_name = 'event.txt'
            f = open(os.path.join(next_dir, file_name), 'w')
            f.write(str(event['gene_node_name']) + ',' + str(event['distance_to_gene_node']) + ',' + str(event['type']) + ',' + str(index))
            f.close()
            gene_subtree.dt_subtree(coalescent_process=species_subtree_coal_process, events=gene_subtree_events, path=next_dir)

        if (event['type'] == 'duplication'):
            Debug.log(header='\n\n\n' + '='*80 + '\nCurrent event:' + '\n',
                      bodies=[event], pformat=True)
            species_skbio_tree = self.species_tree.skbio_tree
            name = self.species_tree.nodes_id_dict[node_id].name

            subtree = species_skbio_tree.find(name).deepcopy()
            subtree_names = [node.name for node in subtree.traverse()]
            subtree_nodes = [node for node in self.species_tree.nodes if node.name in subtree_names]

            species_subtree = SpeciesTree(nodes=subtree_nodes)
            species_subtree.skbio_tree = subtree
            Debug.log(header='\nspecies_subtree_nodes:\n', bodies=species_subtree.nodes)

            distance_above_root = event['distance_to_gene_node'] + coal_distance
            sub_leaves = [int(node_id) for node_id in event['gene_node_name'].strip().split('*')[:-1]]
            Debug.log(header='\nspecies_subtree_coal:\n')
            if (GeneTree.recombination == 1):
                if (GeneTree.hemiplasy == 1):
                    species_subtree_coal_process, chosen_gene_name = species_subtree.incomplete_coalescent(distance_above_root, recombination=1)
                elif (GeneTree.hemiplasy == 0):
                    species_subtree_coal_process = species_subtree.bounded_coalescent(distance_above_root, recombination=1)       
                # original_gene_name = event['gene_node_name']
                # chosen_gene_set = set(chosen_gene_name.split('*')[:-1])
                # original_gene_set = set(original_gene_name.split('*')[:-1])
                # union_gene_set = original_gene_set.union(chosen_gene_set)
                # print(union_gene_set)
            elif (GeneTree.dup_recombination == 0):
                if (GeneTree.hemiplasy == 1):
                    species_subtree_coal_process, chosen_gene_name = species_subtree.incomplete_coalescent(distance_above_root, recombination=0, sub_leaves=sub_leaves)
                elif (GeneTree.hemiplasy == 0):
                    species_subtree_coal_process = species_subtree.bounded_coalescent(distance_above_root, recombination=0, sub_leaves=sub_leaves)        

            Debug.log(header='\nspecies_subtree_coal_process:\n',
                      bodies=[species_subtree_coal_process], pformat=True)

            species_subtree_time_seq = species_subtree.time_sequences(coalescent_process=species_subtree_coal_process)
            Debug.log(header='\nspecies_subtree_time_seq:\n',
                      bodies=[species_subtree_time_seq], pformat=True)

            # save subtree
            # Debug.save_tree_nodes(nodes=species_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'dup', node_id, distance_above_root), 
            #                       distance=distance_above_root)
            
            gene_subtree = GeneTree(time_sequences=species_subtree_time_seq, 
                                    species_tree=species_subtree, 
                                    coalescent_process=species_subtree_coal_process)
            gene_subtree.skbio_tree.length = event['event_height'] - gene_subtree.total_distance
            Debug.log(header='\ngene_subtree nodes:\n', bodies=gene_subtree.nodes)
            # Debug.save_tree_nodes(nodes=gene_subtree.nodes, 
            #                       path=Debug.subtree_file_name('output/subtrees', 'dup', node_id, distance_above_root), 
            #                       mode='a')

            # Debug.save_output(contents=[gene_subtree.skbio_tree],
            #                   path=Debug.subtree_file_name('output/newick_gene_subtrees', 'dup', node_id, distance_above_root))

            Debug.log(header='\ngene_subtree dlt_process:\n')
            gene_subtree_height = gene_subtree.total_distance
            gene_subtree_events = gene_subtree.dlt_process(event=event, distance=event['event_height'] - gene_subtree_height)
            Debug.log(header='\ngene_subtree events:\n',
                      bodies=[gene_subtree_events], pformat=True)

            index = Utility.increment()
            event['index'] = index
            _id = 'dup_subtree_' + str(index)
            next_dir = os.path.join(path, _id)
            os.mkdir(next_dir)
            file_name = 'event.txt'
            f = open(os.path.join(next_dir, file_name), 'w')
            f.write(str(event['gene_node_name']) + ',' + str(event['distance_to_gene_node']) + ',' + str(event['type']) + ',' + str(index))
            f.close()
            gene_subtree.dt_subtree(coalescent_process=species_subtree_coal_process, events=gene_subtree_events, path=next_dir)
        
        return

    # find all the duplication points on the coalescent tree
    # find the corresponding duplicaion subtree
    # do subtree coalescence to obtain the sub_coalescent_tree
    # find all the duplication points on the sub_coalescent_tree
    # recurse
    def dt_subtree(self, coalescent_process, events, path):
        if (path):
            f = open(os.path.join(path, 'gene_tree.txt'), 'w')
            f.write(str(self.skbio_tree))
            f.close()
            f = open(os.path.join(path, 'species_tree.txt'), 'w')
            f.write(str(self.species_tree.skbio_tree))
            f.close()
            self.find_ils(path)
        if (not events):
            return
        GeneTree.full_events += events
        for event in events:
            a = event['gene_node_name'].split('*')[:-1]
            for i in range(len(a)):
                a[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(a[i])
            # Debug.summary(header=event['type'][0] + '\t'
            #      + str(SpeciesTree.global_species_tree.get_fake_id_from_real_id(event['species_node_id'])) + '\t' 
            #      + str(a) + '\n')

            if (event['type'] == 'duplication'):
                # Debug.event_count['D'] += 1
                node_id = None
                coal_distance = None
                if (coalescent_process):        # non-trivial
                    for k, v in coalescent_process.items():
                        for elem in v:
                            if (event['gene_node_name'] in elem['to_set'] and event['gene_node_name'] not in elem['from_set']):
                                node_id = int(k)
                                coal_distance = elem['distance']
                    if (node_id == None):
                        node_id = int(event['gene_node_name'][:-1])
                        coal_distance = 0
                    self.dt_subtree_recurse(event=event, node_id=node_id, coal_distance=coal_distance, path=path)
                else:       # trivial
                    node_id = int(event['gene_node_name'][:-1])
                    coal_distance = 0
                    self.dt_subtree_recurse(event=event, node_id=node_id, coal_distance=coal_distance, path=path)
            elif (event['type'] == 'transfer'):
                # Debug.event_count['T'] += 1
                trans_target_id = event['target_species_id']
                target_height = SpeciesTree.global_species_tree.distance_to_leaf(trans_target_id, 0)
                distance_above_target = event['event_height'] - target_height
                self.dt_subtree_recurse(event=event, node_id=trans_target_id, coal_distance=distance_above_target, path=path)
            elif (event['type'] == 'loss'):
                # Debug.summary(header=event['type'][0] + '\t'
                #  + str(SpeciesTree.global_species_tree.get_fake_id_from_real_id(event['species_node_id'])) + '\t' 
                #  + str(a) + '\n')
                index = Utility.increment()
                event['index'] = index
                Debug.event_count['L'] += 1
                file_name = 'loss_' + str(event['distance_to_gene_node'])
                f = open(os.path.join(path, file_name), 'w')
                f.write(str(event['gene_node_name']) + ',' + str(event['distance_to_gene_node']) + ',' + str(index))
                f.close()
        return