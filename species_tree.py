from shared_utility import *

import numpy as np
from statistics import mean
from collections import defaultdict


class SpeciesTree(GenericTree):
    # static properties
    global_species_tree = None
    lambda_coal = None

    def __init__(self,
                 newick_path=None,
                 nodes=None):
        GenericTree.__init__(self)
        if (not newick_path):
            self.nodes = nodes
        else:
            self.construct_species_nodes(newick_path)

        max_node_id = -1
        for node in self.nodes:
            if (node.node_id > max_node_id):
                max_node_id = node.node_id
                self.root = node
            self.nodes_id_dict[node.node_id] = node
            self.nodes_name_dict[node.name] = node

            node.clade = []
            for i in range(len(node.name)):  
                char = node.name[i]
                node_id = self.node_by_name(char).node_id
                node.clade.append(node_id)

            node.clade_split = []
            if (node.children and not node.clade_split):
                for i in range(len(node.children)):
                    node_name = self.node_by_id(node.children[i]).name
                    split = []
                    for j in range(len(node_name)):  
                        char = node_name[j]
                        node_id = self.node_by_name(char).node_id
                        split.append(node_id)
                    node.clade_split.append(split)

        self.leaves = [node.node_id for node in self.nodes if not node.children]
        self.total_distance = self.distance_to_root_recurse(node_id=self.leaves[0])
        return

    def construct_species_nodes(self, newick_path):
        output_path = 'output/temp_species_nodes_table.txt'
        self.skbio_tree = super().newick_to_table(input_path=newick_path, output_path=output_path)
        super().construct_nodes(path=output_path, process_tree=True)
        return

    # def assign_fake_id(self):
    #     curr_index = 0
    #     for n in SpeciesTree.global_species_tree.skbio_tree.postorder():
    #         for c in n.children:
    #             node = self.node_by_name(c.name)
    #             node.fake_node_id = curr_index
    #             curr_index += 1

    #     root = self.node_by_name(SpeciesTree.global_species_tree.skbio_tree.name)
    #     root.fake_node_id = curr_index


    def get_lambda_coal(self, clade_set):
        indices = []
        for clade in clade_set:
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
        return mean(SpeciesTree.lambda_coal[indices])

    # checking whether a given clade is in the target set
    # modified for the "*" representation
    def star_in_set(self, target, clade):
        if (len(target) <= len(clade)):
            splited_target = target.split('*')[:-1]
            splited_clade = clade.split('*')[:-1]
            return set(splited_target).issubset(set(splited_clade))
        else:
            return False

    def star_sorted(self, couple):
        string = ''
        for e in couple:
            string += e
        splited = string.split('*')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '*' for e in splited]

    # This is the recursive part of the multi-species coalescent process:
    # Given a set of n genes gathering into a branch in the species tree from the bottom,
    #   whenever we come across a point of coalescence, we randomly merge 2 elements in the gene sets,
    #   and record the set before the new coalescence, named "from_set", and the set after the coalescence,
    #   named "to_set", and the distance from the last coalescent event or the bottom of the branch.
    def coalescent_recurse(self, node_id, distance, clade_set, coalescent_process):
        if (len(clade_set[node_id]) <= 1):
            return clade_set[node_id]
        else:
            lambda_c = len(clade_set[node_id]) * self.get_lambda_coal(clade_set[node_id])    # rate of coalescence
            distance_fake = np.random.exponential(scale=1.0/lambda_c)
            if (distance < distance_fake):      # no coalescent event anymore in this branch
                return clade_set[node_id]
            else:
                if (len(clade_set[node_id]) >= 2):   # when coalescent, randomly merge 2 elements in the gene sets
                    temp_set = sorted(clade_set[node_id])
                    couple = np.random.choice(clade_set[node_id], size=2, replace=False)
                    clade_set[node_id] = [''.join(self.star_sorted(couple))] + [e for e in clade_set[node_id] if e not in couple]

                    # print process
                    Debug.log(header="initial node " + str(node_id) + ": " + str(temp_set) + '\n')
                    Debug.log(header="coalescent at node " + str(node_id) + ": " + str(clade_set[node_id]) + ", " + "distance = " + str(distance_fake) + '\n')

                    # save process
                    coalescent_process[str(node_id)].append({
                        'from_set': temp_set, 
                        'to_set': clade_set[node_id].copy(),
                        'distance': distance_fake
                    })
                else:
                    return clade_set[node_id]     # stop when gene set only has one single element
                distance = distance - distance_fake
                self.coalescent_recurse(node_id=node_id, 
                                        distance=distance, 
                                        clade_set=clade_set,
                                        coalescent_process=coalescent_process)     # use recursion to simulate the case when there is more than one coalescent events in the branch
        return clade_set[node_id]

    # the main multi-species coalecent function
    def coalescent(self, distance_above_root):
        nodes = self.nodes
        root = self.root
        coalescent_process = defaultdict(list)
        clade_set_into_root = None

        old_leaves = [node.node_id for node in nodes if not node.children]      # leaves of the given species tree
        new_leaves = []     # leaves set will be updated in the loop
        clade_set = {}      # set of extant species that an ancestral gene will finally be fixed in
        labelled = {}       # avoid doing repeated coalescence
        for node in nodes:
            labelled[node.node_id] = False
            clade_set[node.node_id] = [str(node.node_id) + '*'] if not node.children else []

        while (True):
            for leaf in old_leaves:
                if (leaf == root.node_id):
                    clade_set_into_root = self.coalescent_recurse(node_id=root.node_id, 
                                                                  distance=distance_above_root, 
                                                                  clade_set=clade_set,
                                                                  coalescent_process=coalescent_process)
                    break
                else:
                    parent = self.nodes_id_dict[leaf].parent_id
                    children = self.nodes_id_dict[parent].children
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(clade_set[children[0]]) != 0 
                        and len(clade_set[children[1]]) != 0):
                        self.coalescent_recurse(node_id=children[0], 
                                                distance=self.nodes_id_dict[children[0]].distance_to_parent,
                                                clade_set=clade_set,
                                                coalescent_process=coalescent_process)
                        self.coalescent_recurse(node_id=children[1], 
                                                distance=self.nodes_id_dict[children[1]].distance_to_parent,
                                                clade_set=clade_set,
                                                coalescent_process=coalescent_process)
                        # the clade set of the parent before coalescence is the union of the clade set of its children after coalescence                        
                        clade_set[parent] = list(set().union(clade_set[children[0]], clade_set[children[1]]))    
                        if (len(new_leaves) > 0):
                            new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                        new_leaves.append(parent)
                    else:
                        new_leaves.append(leaf)     # updating leaves set
            if (leaf == root.node_id):
                break
            temp_new_leaves = []
            for new_leaf in new_leaves:
                if (new_leaf not in temp_new_leaves):
                    temp_new_leaves.append(new_leaf)
            old_leaves = temp_new_leaves.copy()
            new_leaves = []
            labelled = {}
            for node in nodes:
                labelled[node.node_id] = False
        return coalescent_process, clade_set_into_root

    # this is the main function of doing sub_species-tree coalescence
    # this will be used when modelling duplications and transfers
    # Given a sub_species_tree, and a subset of leaves in the sub_species_tree, named "sub_leaves"
    #   this function dose a multi-species coalescence based on the sub_species-tree only considering the sub_leaves
    def sub_leaves_coalescent(self, distance_above_root, sub_leaves):
        nodes = self.nodes
        root = self.root
        coalescent_process = defaultdict(list)

        old_leaves = [node.node_id for node in nodes if not node.children]
        new_leaves = []
        clade_set = {}
        labelled = {}
        mark = {} 

        for node in nodes:
            mark[node.node_id] = True
            labelled[node.node_id] = False
            clade_set[node.node_id] = [str(node.node_id) + '*'] if not node.children else [] 
        for leaf in old_leaves: mark[leaf] = False
        for leaf in sub_leaves: mark[leaf] = True       # mark the leaves of the sub_species_tree but not in the sub_leaves set as FALSE

        while (True):
            for leaf in old_leaves:
                if (leaf == root.node_id):
                    clade_set_into_root = self.coalescent_recurse(node_id=root.node_id, 
                                            distance=distance_above_root, 
                                            clade_set=clade_set,
                                            coalescent_process=coalescent_process)
                    break
                else:
                    parent = self.nodes_id_dict[leaf].parent_id
                    children = self.nodes_id_dict[parent].children
                    if (labelled[leaf]):
                        continue
                    labelled[leaf] = True
                    if (len(clade_set[children[0]]) != 0 
                        and len(clade_set[children[1]]) != 0):
                            if (mark[children[0]] and mark[children[1]]):       # when both children are marked, do normal coalescence
                                self.coalescent_recurse(node_id=children[0], 
                                                        distance=self.nodes_id_dict[children[0]].distance_to_parent,
                                                        clade_set=clade_set,
                                                        coalescent_process=coalescent_process)
                                self.coalescent_recurse(node_id=children[1], 
                                                        distance=self.nodes_id_dict[children[1]].distance_to_parent,
                                                        clade_set=clade_set,
                                                        coalescent_process=coalescent_process)
                                clade_set[parent] = list(set().union(clade_set[children[0]], clade_set[children[1]]))
                                if (len(new_leaves) > 0):
                                    new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                                new_leaves.append(parent)
                            elif not (mark[children[0]] or mark[children[1]]):      # when neither childer is marked, we mark the parent as FALSE
                                mark[parent] = False
                                clade_set[parent] = '$'
                                if (len(new_leaves) > 0):
                                    new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                                new_leaves.append(parent)
                            elif (mark[children[0]]):       # when only one child is marked, let the clade set of the unmarked child be empty
                                clade_set[parent] = clade_set[children[0]]
                                if (len(new_leaves) > 0):
                                    new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                                new_leaves.append(parent)
                            elif (mark[children[1]]):
                                clade_set[parent] = clade_set[children[1]]
                                if (len(new_leaves) > 0):
                                    new_leaves = [e for e in new_leaves if e != children[0] and e != children[1]]
                                new_leaves.append(parent)
                    else:
                        new_leaves.append(leaf)
            if (leaf == root.node_id):
                break
            temp_new_leaves = []
            for new_leaf in new_leaves:
                if (new_leaf not in temp_new_leaves):
                    temp_new_leaves.append(new_leaf)
            old_leaves = temp_new_leaves.copy()
            new_leaves = []
            labelled = {}
            for node in nodes:
                labelled[node.node_id] = False
        return coalescent_process, clade_set_into_root

    def incomplete_coalescent(self, distance_above_root, recombination=1, sub_leaves=None): # hemiplasy
        if (recombination == 1):
            full_coal_process, genes_into_root = self.coalescent(distance_above_root)
        elif (recombination == 0):
            full_coal_process, genes_into_root = self.sub_leaves_coalescent(distance_above_root, sub_leaves)
        chosen_gene = np.random.choice(genes_into_root)
        sub_coal_process = self.filter_coal_process(full_coal_process, chosen_gene)
        return sub_coal_process, chosen_gene

    def bounded_coalescent(self, distance_above_root, recombination=1, sub_leaves=None): # no hemiplasy
        if (recombination == 0):
            coal_process, genes_into_root = self.sub_leaves_coalescent(distance_above_root, sub_leaves)
        elif (recombination == 1):
            coal_process, genes_into_root = self.coalescent(distance_above_root)
        if (len(genes_into_root) == 1):
            return coal_process
        else:
            return self.bounded_coalescent(distance_above_root, recombination, sub_leaves)

    # given a coalescent process obtained by incomplete coalescent,
    # one may have more than one subtrees in the full_coal_process,
    # we can choose a subtree rooted at the chosen_gene,
    # at this stage, the subtree is represented as a modified coal_process.
    def filter_coal_process(self, full_coal_process, chosen_gene):
        coal_process = defaultdict(list)
        for k, v in full_coal_process.items():
            for elem in v:
                distance = elem['distance']
                from_set = []
                to_set = []
                for clade in elem['from_set']:
                    if (self.star_in_set(target=clade, clade=chosen_gene)):
                        from_set.append(clade)
                for clade in elem['to_set']:
                    if (self.star_in_set(target=clade, clade=chosen_gene)):
                        to_set.append(clade)
                if (to_set):
                    coal_process[k].append({
                        'from_set': from_set, 
                        'to_set': to_set,
                        'distance': distance
                    })
        return coal_process
    
    # find the ancestors of the given leaf in reverse time order
    def find_ancestors(self, leaf_name, coalescent_process):
        sequence = []
        for k, v in coalescent_process.items():
            branch_distance = 0.0
            for elem in v:
                branch_distance += elem['distance']
                if (leaf_name in elem['from_set'] and leaf_name not in elem['to_set']):
                    for e in elem['to_set']:
                        if len(leaf_name) < len(e) and self.star_in_set(leaf_name, e):
                            # species_node_id = int(k)
                            # species_node_height = super().distance_to_leaf(int(k), branch_distance=0)
                            coal_height = super().distance_to_leaf(int(k), branch_distance=branch_distance)
                            # pair = (ancestor, coal_height)
                            pair = (e, coal_height)
                            sequence.append(pair)
                            sequence += self.find_ancestors(leaf_name=e, 
                                                            coalescent_process=coalescent_process)
        return sequence

    # backward-in-time coalescent process
    # modified data structure for constructing the coalescent tree in newick format
    def time_sequences(self, coalescent_process):
        time_sequences = {}
        for leaf in self.leaves:
            time_sequences[str(leaf)] = self.find_ancestors(leaf_name=str(leaf)+'*', 
                                                            coalescent_process=coalescent_process)
        return time_sequences