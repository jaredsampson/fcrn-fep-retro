'''
My awesome script.
'''
import argparse
import logging
import random
from typing import Optional

import networkx as nx

from schrodinger.application.desmond.constants import WT_IDENTIFIER
from schrodinger.application.desmond.file_utils import get_graph_file_path_base
from schrodinger.application.scisol.fep import graph
# from schrodinger.application.scisol.fep.cycle_closure import CycleFinder
from schrodinger.utils import log


# Configure logging
logger = log.get_output_logger(__name__)
logger.level = logging.WARNING


DEFAULT_SEED = 0
DEFAULT_MAX_PATHS = 2
MAX_CYCLE_SIZE = 4
MIN_EDGES_PER_NODE = 2


# def count_cycles_with_node(ggraph: graph.Graph, node: graph.Node,
#                            max_visits=100000) -> collections.Counter:
#     """
#     Count number of cycles in `ggraph` containing `node`, grouped by length.

#     This routine counts the lengths of (possibly a subset of) cycles of a
#     given graph.

#     For example, if a graph has 2 cycles of length 3. This function will
#     return a data structure similar to {3: 2}.

#     :param graph: Graph to inspect.
#     :param node: Node to check for cycles.
#     :param max_visits: Threshold on the amount of work that can be done to
#         obtain cycles. By default, we evaluate one order of magnitude less than
#         the default for CycleFinder to increase speed.
#     :return: A dictionary mapping cycle lengths to number
#         of occurrences of cycles of thos lengths
#     """

#     def canonicalize(cycle: list[graph.Node]) -> tuple[str, ...]:
#         return tuple(sorted(node.id for node in cycle))

#     cycle_finder = CycleFinder(ggraph, max_visits=max_visits)

#     # Make sure we're not double-counting permutations of the same cycle.
#     cycles = collections.defaultdict(set)
#     for cycle in cycle_finder.find_cycles():
#         canonical_cycle = canonicalize(cycle)
#         if node.id in canonical_cycle:
#             cycles[len(canonical_cycle)].add(canonical_cycle)

#     # Now count the number of cycles grouping by length.
#     return collections.Counter({length: len(v) for length, v in cycles.items()})


def count_simple_paths(ggraph, n0, n1, depth):
    '''Count the number of simple paths between two nodes.'''
    paths = nx.all_simple_paths(ggraph, n0, n1, cutoff=depth)
    return len(list(paths))


def make_sparse_graph(ggraph: graph.Graph, max_paths: int, max_cycle_size: int,
                      seed: Optional[int] = None) -> graph.Graph:
    '''Remove redundant edges to make a sparse graph.'''
    logger.info('\n- Removing redundant edges...')

    removed_edge_count = 0
    edges = list(ggraph.edges())
    # Shuffle edges to avoid bias in the order of removal
    random.seed(seed)
    for edge in random.sample(edges, k=len(edges)):
        n0, n1 = edge

        # An edge is considered redundant if it connects two nodes that are
        # already connected by more than `max_paths` distinct simple paths less
        # than `max_cycle_size`, and both nodes are overconnected, meaning they
        # each have more than `MIN_EDGES_PER_NODE` edges.
        depth = max_cycle_size - 1
        too_many_paths = count_simple_paths(ggraph, n0, n1, depth) > max_paths
        n0_overconnected = len(n0.edges()) > MIN_EDGES_PER_NODE
        n1_overconnected = len(n1.edges()) > MIN_EDGES_PER_NODE
        is_redundant = too_many_paths and n0_overconnected and n1_overconnected

        # Keep all edges that are direct mutations from WT node
        is_direct_edge = WT_IDENTIFIER in [n0.struc.title, n1.struc.title]

        if is_redundant and not is_direct_edge:
            title = edge.short_id_title
            ggraph.remove_edge(edge)
            logger.info(f'  - removed edge: {title}')
            removed_edge_count += 1
        else:
            logger.debug(f'  - kept edge: {edge.short_id_title}')

    logger.warning(f'- Removed {removed_edge_count} redundant edges')

    return ggraph


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('fmp',
                        help='path to FMP file')
    parser.add_argument('-max-paths', type=int, default=DEFAULT_MAX_PATHS,
                        help='maximum number of paths to keep between nodes')
    parser.add_argument('-max-cycle-size', type=int, default=MAX_CYCLE_SIZE,
                        help='maximum size of cycles to consider; larger '
                        'values allow more edges to be removed')
    parser.add_argument('-seed', type=int, default=DEFAULT_SEED,
                        help='random seed for reproducibility')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('-debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


def main(argv=None):
    '''Main workflow, optionally callable like subprocess with list of args.'''
    args = parse_args(argv)

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    logger.info(f'\n- Loading graph from {args.fmp}')
    ggraph = graph.Graph.deserialize(args.fmp)
    logger.info(f'  - Input graph has {ggraph.number_of_nodes()} nodes and '
                f'{ggraph.number_of_edges()} edges')

    sparse_ggraph = make_sparse_graph(ggraph,
                                      max_paths=args.max_paths,
                                      max_cycle_size=args.max_cycle_size,
                                      seed=args.seed)

    fmp_base = get_graph_file_path_base(args.fmp)
    out_fmp = f'{fmp_base}_sparse.fmp'
    sparse_ggraph.write(out_fmp)
    logger.warning(f'- Wrote {out_fmp}')
    logger.info(f'  - Output graph has {sparse_ggraph.number_of_nodes()} nodes '
                f'and {sparse_ggraph.number_of_edges()} edges')


if __name__ == '__main__':
    main()
