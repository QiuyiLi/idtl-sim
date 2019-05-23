import skbio
import pprint


class Utility(object):
    count = 0

    @staticmethod
    def increment():
        Utility.count += 1
        return Utility.count


class Debug(object):
    log_file = None
    summary_file = None
    recon_file = None
    event_count = {'D': 0, 'L': 0, 'T': 0, 'I': 0, 'S':0, 'SL': 0, 'TL': 0, 'IL': 0, 'DL': 0}

    def __init__(self):
        pass

    @staticmethod
    def subtree_file_name(path, event, node_id, distance):
        return '{}/subtree_{}_{:03d}_{:07d}.txt'.format(path, event, node_id, int(distance * 10000))

    @staticmethod
    def save_tree_nodes(nodes, path, mode='w', distance=None):
        f = open(path, mode)
        if (distance): f.write(str(distance) + '\n')
        f.write('id\tname\tparent\td2p\n')
        for node in nodes:
            f.write(str(node.node_id) + '\t' + node.name + '\t' + str(node.parent_id) + '\t' + str(node.distance_to_parent) + '\n')
        f.close()

    @staticmethod
    def save_output(contents, path, mode='w'):
        f = open(path, mode)
        for content in contents:
            f.write(str(content))
            f.write('\n')
        f.close()

    @staticmethod
    def log(header, bodies=[], pformat=False):
        Debug.log_file.write(header)
        for body in bodies:
            if (pformat):
                Debug.log_file.write(pprint.pformat(body))
            else:
                Debug.log_file.write(str(body))
            Debug.log_file.write('\n')

    @staticmethod
    def summary(header, bodies=[], pformat=False):
        Debug.summary_file.write(header)
        for body in bodies:
            if (pformat):
                Debug.summary_file.write(pprint.pformat(body))
            else:
                Debug.summary_file.write(str(body))
            Debug.summary_file.write('\n')
    
    @staticmethod
    def recon_table(header, bodies=[], pformat=False):
        Debug.recon_file.write(header)
        for body in bodies:
            if (pformat):
                Debug.recon_file.write(pprint.pformat(body))
            else:
                Debug.recon_file.write(str(body))
            Debug.recon_file.write('\n')


class TreeNode(object):
    def __init__(self,
                 node_id=None,
                 name=None,
                 parent=None,
                 distance_to_parent=None,
                 children=None):
        self.node_id = node_id
        self.fake_node_id = -1
        self.name = name
        self.parent_id = parent
        self.distance_to_parent = distance_to_parent
        self.children = children if children else []
        self.distance_to_children = []
        self.clade = []
        self.clade_split = []
        self.events = []
        return
    
    def __repr__(self):
        return "node_id: {}, name: {}, parent_id: {}, distance_to_parent: {}, children: {}, distance_to_children: {}, clade: {}, clade_split: {}".format(
                self.node_id, 
                self.name,
                self.parent_id, 
                self.distance_to_parent,
                self.children,
                self.distance_to_children,
                self.clade,
                self.clade_split)


class GenericTree(object):
    def __init__(self):
        self.skbio_tree = None
        self.nodes = []
        self.root = None
        self.leaves = []
        self.nodes_id_dict = {}
        self.nodes_name_dict = {}
        self.total_distance = -1
        return

    def newick_to_table(self, output_path, input_path=None, skbio_tree=None):
        tree = None
        if (skbio_tree):
            tree = skbio_tree
        else:
            f = open(input_path)
            tree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
            f.close()

        def _parse(tree):
            node = {
                'object': tree,
                'name': tree.name,
                'parent': tree.parent,
                'children': [],
                'distance': tree.length
            }
            if tree.is_tip():
                return node
            for children in tree.children:
                node['children'].append(_parse(children))
            return node

        def _rename(node):
            ret = node.copy()
            name = ''
            for i in range(len(ret['children'])):
                ret['children'][i] = _rename(ret['children'][i])
                name += ret['children'][i]['name']
            if (not ret['name']):
                ret['name'] = name
            return ret

        def _to_list(node, root):
            d = node.copy()
            ret = []
            for i in range(len(d['children'])):
                ret += _to_list(d['children'][i], root=root)
            del d['children']
            d['distance_to_root'] = d['object'].distance(root)
            ret.append(d)
            if (d['object'] is root):
                ret = sorted(ret, key=lambda x: x['distance_to_root'], reverse=True)
            return ret

        root = _parse(tree)
        root = _rename(root)
        nodes = _to_list(root, root=root['object'])

        # output to file
        f = open(output_path, 'w')
        f.write('id\tname\tparent\td2p\n')
        for i in range(len(nodes)):
            parent_index = 'None'
            for j in range(len(nodes)):
                if (nodes[i]['parent'] is nodes[j]['object']):
                    parent_index = str(j)
            f.write(str(i) + '\t' + nodes[i]['name'] + '\t' + parent_index + '\t' + str(nodes[i]['distance']) + '\n')
        f.close()

        return tree

    def process_tree_recurse(self, tree):
        if (tree.name):
            return tree.name
        else:
            name = ''
            for child in tree.children:
                name += self.process_tree_recurse(child)
            tree.name = name
            return tree.name
    
    def construct_nodes(self, path, process_tree=False):
        # open table file and construct tree nodes
        f = open(path, 'r')
        f.readline()
        for line in f:
            splited = line.strip().split('\t')
            node_id = int(splited[0])
            name = splited[1]
            parent = splited[2]
            if (parent == 'None'):
                parent = -1
            d2p = splited[3]
            if (d2p == 'None'):
                d2p = -1.0
            self.nodes.append(TreeNode(node_id=node_id,
                                       name=name,
                                       parent=int(parent),
                                       distance_to_parent=float(d2p)))
        # process tree
        if (process_tree):
            self.process_tree_recurse(self.skbio_tree)
        # create dict
        for node in self.nodes:
            self.nodes_name_dict[node.name] = node
        # find children
        for node in self.nodes:
            children, node.distance_to_children = self.children_distances(self.skbio_tree, node.name)
            for child in children:
                node_id = self.nodes_name_dict[child.name].node_id
                node.children.append(node_id)
        return

    def children_distances_recurse(self, tree, name):
        ret = None
        if (tree.name == name):
            distances = []
            for child in tree.children:
                distances.append(tree.distance(child))
            ret = (tree.children, distances)
            return ret
        else:
            for child in tree.children:
                ret = self.children_distances_recurse(child, name)
                if (ret): return ret

    def children_distances(self, tree, name):
        return self.children_distances_recurse(tree, name)

    def node_by_id(self, node_id):
        return self.nodes_id_dict[node_id]

    def node_by_name(self, name):
        return self.nodes_name_dict[name]

    # def fake(self):
    #     temp_nodes = self.nodes.copy()
    #     curr_index = 0
    #     while (len(temp_nodes) > 0):
    #         parent_nodes = []
    #         for node in temp_nodes:
    #             if (node.fake_node_id == -1):
    #                 node.fake_node_id = curr_index
    #                 curr_index += 1
    #                 if (node.parent_id != -1):
    #                     parent_nodes.append(self.node_by_id(node.parent_id))
    #         temp_nodes = parent_nodes

    def post_order_fake_id_recurse(self, index, node):
        curr_index = index
        for child in node.children:
            child_node = self.node_by_id(child)
            curr_index = self.post_order_fake_id_recurse(curr_index, child_node)
        node.fake_node_id = curr_index
        curr_index += 1
        return curr_index

    def post_order_fake_id(self):
        root = self.root
        self.post_order_fake_id_recurse(0, root)

    def get_fake_id_from_real_id(self, real_id):
        node = self.node_by_id(int(real_id))
        fake_id = node.fake_node_id
        return fake_id
    
    # find the distance of a given node to the root
    # needed when finding the walking distance
    def distance_to_root_recurse(self, node_id):        
        if (self.nodes_id_dict[node_id].parent_id < 0 or 
            node_id == self.root.node_id):
            return 0
        else:
            d2p = self.nodes_id_dict[node_id].distance_to_parent
            parent = self.nodes_id_dict[node_id].parent_id
            return d2p + self.distance_to_root_recurse(parent)

    # given a coalescent event happening at "branch_distance" above a speices node with "node_id"
    # find the distance of this event to the bottom of the tree
    # needed when assigning ids to the coalescent tree
    def distance_to_leaf(self, node_id, branch_distance):
        return branch_distance + (self.total_distance - self.distance_to_root_recurse(node_id))